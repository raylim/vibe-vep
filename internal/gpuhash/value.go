// Package gpuhash provides a GPU-resident (or in-memory) open-addressing hash
// table for batch annotation source lookups (AlphaMissense, ClinVar, SIGNAL).
//
// The hot path in the annotation pipeline enriches one variant at a time with
// sequential SQLite prepared-statement executions, which are ~1–5 µs each and
// run single-threaded in the ordered collector callback.  This package replaces
// that with a single batch lookup call that is either:
//
//   - An in-memory Go map (non-cuda build): ~150–300 ns/lookup, fully parallel
//   - A CUDA kernel on the GPU hash table (cuda build): ~50 ms per 1 M variants
//
// # Key encoding
//
// Lookups use a 64-bit FNV-1a hash of the canonical key
// "chrom\x00pos(le32)\x00ref\x00alt".  Hash collisions in a 10 M-entry table
// are vanishingly rare (~2×10⁻⁷ probability per variant); the CPU path does
// full-key equality checks, while the GPU path accepts the negligible error.
//
// # Value encoding
//
// Values are stored in compact fixed-width structs to maximise GPU cache
// efficiency.  Variable-length text fields (ClinVar disease name, SIGNAL freq
// string) are decoded back to strings by the caller after lookup.  The Slot
// size is 32 bytes, giving natural alignment on all 32-bit warp lanes.
package gpuhash

// AMClass encodes the AlphaMissense pathogenicity class.
type AMClass uint8

const (
	AMClassAbsent           AMClass = 0
	AMClassLikelyBenign     AMClass = 1
	AMClassAmbiguous        AMClass = 2
	AMClassLikelyPathogenic AMClass = 3
)

// AMClassString returns the canonical string for an AlphaMissense class.
var AMClassString = [4]string{"", "likely_benign", "ambiguous", "likely_pathogenic"}

// EncodeAMClass converts an AM class string to its compact encoding.
func EncodeAMClass(s string) AMClass {
	switch s {
	case "likely_benign":
		return AMClassLikelyBenign
	case "ambiguous":
		return AMClassAmbiguous
	case "likely_pathogenic":
		return AMClassLikelyPathogenic
	}
	return AMClassAbsent
}

// CVSig encodes ClinVar clinical significance as a uint8.
type CVSig uint8

const (
	CVSigAbsent                 CVSig = 0
	CVSigPathogenic             CVSig = 1
	CVSigLikelyPathogenic       CVSig = 2
	CVSigPathogenicOrLikely     CVSig = 3
	CVSigBenign                 CVSig = 4
	CVSigLikelyBenign           CVSig = 5
	CVSigBenignOrLikely         CVSig = 6
	CVSigUncertainSignificance  CVSig = 7
	CVSigConflicting            CVSig = 8
	CVSigNotProvided            CVSig = 9
	CVSigDrugResponse           CVSig = 10
	CVSigRiskFactor             CVSig = 11
	CVSigOther                  CVSig = 255
)

// CVSigStrings maps CVSig codes to their canonical string representation.
var CVSigStrings = map[CVSig]string{
	CVSigPathogenic:            "Pathogenic",
	CVSigLikelyPathogenic:      "Likely_pathogenic",
	CVSigPathogenicOrLikely:    "Pathogenic/Likely_pathogenic",
	CVSigBenign:                "Benign",
	CVSigLikelyBenign:          "Likely_benign",
	CVSigBenignOrLikely:        "Benign/Likely_benign",
	CVSigUncertainSignificance: "Uncertain_significance",
	CVSigConflicting:           "Conflicting_interpretations_of_pathogenicity",
	CVSigNotProvided:           "not_provided",
	CVSigDrugResponse:          "drug_response",
	CVSigRiskFactor:            "risk_factor",
}

// EncodeCVSig converts a ClinVar CLNSIG string to its compact encoding.
func EncodeCVSig(s string) CVSig {
	switch s {
	case "Pathogenic":
		return CVSigPathogenic
	case "Likely_pathogenic":
		return CVSigLikelyPathogenic
	case "Pathogenic/Likely_pathogenic", "Pathogenic,_Likely_pathogenic":
		return CVSigPathogenicOrLikely
	case "Benign":
		return CVSigBenign
	case "Likely_benign":
		return CVSigLikelyBenign
	case "Benign/Likely_benign", "Benign,_Likely_benign":
		return CVSigBenignOrLikely
	case "Uncertain_significance", "Uncertain significance":
		return CVSigUncertainSignificance
	case "Conflicting_interpretations_of_pathogenicity",
		"Conflicting interpretations of pathogenicity":
		return CVSigConflicting
	case "not_provided", "Not provided":
		return CVSigNotProvided
	case "drug_response", "Drug response":
		return CVSigDrugResponse
	case "risk_factor", "Risk factor":
		return CVSigRiskFactor
	case "":
		return CVSigAbsent
	}
	return CVSigOther
}

// CVRevStat encodes ClinVar review status.
type CVRevStat uint8

const (
	CVRevStatAbsent                  CVRevStat = 0
	CVRevStatNoCriteria              CVRevStat = 1
	CVRevStatSingleSubmitter         CVRevStat = 2
	CVRevStatMultipleSubmitters      CVRevStat = 3
	CVRevStatReviewedByExpertPanel   CVRevStat = 4
	CVRevStatPracticeGuideline       CVRevStat = 5
	CVRevStatOther                   CVRevStat = 255
)

// CVRevStatStrings maps CVRevStat codes to canonical strings.
var CVRevStatStrings = map[CVRevStat]string{
	CVRevStatNoCriteria:            "no_assertion_criteria_provided",
	CVRevStatSingleSubmitter:       "criteria_provided,_single_submitter",
	CVRevStatMultipleSubmitters:    "criteria_provided,_multiple_submitters,_no_conflicts",
	CVRevStatReviewedByExpertPanel: "reviewed_by_expert_panel",
	CVRevStatPracticeGuideline:     "practice_guideline",
}

// EncodeCVRevStat converts a ClinVar review status string to its encoding.
func EncodeCVRevStat(s string) CVRevStat {
	switch s {
	case "no_assertion_criteria_provided":
		return CVRevStatNoCriteria
	case "criteria_provided,_single_submitter":
		return CVRevStatSingleSubmitter
	case "criteria_provided,_multiple_submitters,_no_conflicts":
		return CVRevStatMultipleSubmitters
	case "reviewed_by_expert_panel":
		return CVRevStatReviewedByExpertPanel
	case "practice_guideline":
		return CVRevStatPracticeGuideline
	case "":
		return CVRevStatAbsent
	}
	return CVRevStatOther
}

// SigStatus encodes SIGNAL mutation status.
type SigStatus uint8

const (
	SigStatusAbsent   SigStatus = 0
	SigStatusGermline SigStatus = 1
	SigStatusSomatic  SigStatus = 2
)

// SigStatusStrings maps SigStatus to canonical strings.
var SigStatusStrings = [3]string{"", "germline", "somatic"}

// EncodeSigStatus converts a SIGNAL mutation status string.
func EncodeSigStatus(s string) SigStatus {
	switch s {
	case "germline", "Germline":
		return SigStatusGermline
	case "somatic", "Somatic":
		return SigStatusSomatic
	}
	return SigStatusAbsent
}

// Value is the decoded result of a hash table lookup.
// All fields are zero-valued when absent.
type Value struct {
	AMScore    float32
	AMClass    AMClass
	CVSig      CVSig
	CVRevStat  CVRevStat
	CVClnDN    string // disease name; retrieved separately (variable-length)
	SigStatus  SigStatus
	SigCount   uint32
	SigFreq    float32
}

// HasAM returns true if an AlphaMissense score is present.
func (v Value) HasAM() bool { return v.AMClass != AMClassAbsent }

// HasCV returns true if a ClinVar entry is present.
func (v Value) HasCV() bool { return v.CVSig != CVSigAbsent }

// HasSig returns true if a SIGNAL entry is present.
func (v Value) HasSig() bool { return v.SigStatus != SigStatusAbsent }

// Slot is the fixed-width hash table entry stored in both GPU memory and the
// CPU in-memory table.  The layout is packed to 32 bytes for CUDA alignment.
//
//	Offset  Size  Field
//	 0       8    key_hash (uint64; 0 = empty slot)
//	 8       4    am_score (float32)
//	12       1    am_class (AMClass uint8)
//	13       1    cv_sig   (CVSig uint8)
//	14       1    cv_revstat (CVRevStat uint8)
//	15       1    sig_status (SigStatus uint8)
//	16       4    sig_count (uint32)
//	20       4    sig_freq (float32)
//	24       8    reserved / padding
type Slot struct {
	Hash      uint64
	AMScore   float32
	AMClass   AMClass
	CVSig     CVSig
	CVRevStat CVRevStat
	SigStatus SigStatus
	SigCount  uint32
	SigFreq   float32
	_         [8]byte // padding to 32 bytes
}
