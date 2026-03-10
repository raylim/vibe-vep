package clinvar

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"regexp"
	"strconv"
	"strings"
)

// SummaryEntry holds a single ClinVar variant_summary.txt.gz record
// filtered and normalized for benchmarking.
type SummaryEntry struct {
	Chrom       string // e.g. "17"
	Pos         int64  // VCF 1-based position
	Ref         string // reference allele
	Alt         string // alternate allele
	Gene        string // gene symbol
	Transcript  string // RefSeq NM_ transcript (unversioned, e.g. "NM_000546")
	Protein     string // protein notation e.g. "p.Arg273Cys"
	ClnSig      string // e.g. "Pathogenic"
	RevStatus   string // e.g. "reviewed by expert panel"
	VariantType string // "snv", "deletion", "insertion", "indel"
	IsMANE      bool   // set later by caller after MANE lookup
}

var reProtein = regexp.MustCompile(`\(p\.([^)]+)\)`)
var reTranscript = regexp.MustCompile(`^([A-Z]{2}_\d+)`)

// ParseSummaryFile reads a (optionally gzipped) ClinVar variant_summary.txt[.gz]
// and returns all GRCh38 pathogenic/likely-pathogenic SNVs and small indels
// (deletions, insertions, indels) that carry a protein-level HGVS annotation.
func ParseSummaryFile(path string) ([]SummaryEntry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open clinvar summary: %w", err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("gzip clinvar summary: %w", err)
		}
		defer gz.Close()
		r = gz
	}

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 4*1024*1024), 4*1024*1024)

	if !scanner.Scan() {
		return nil, fmt.Errorf("clinvar summary: empty file")
	}
	cols := indexCols(strings.Split(scanner.Text(), "\t"),
		"Type", "Name", "GeneSymbol", "ClinicalSignificance", "ReviewStatus",
		"Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")

	required := []string{"Type", "Name", "GeneSymbol", "ClinicalSignificance",
		"Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"}
	for _, c := range required {
		if cols[c] < 0 {
			return nil, fmt.Errorf("clinvar summary: missing column %q", c)
		}
	}

	var entries []SummaryEntry
	for scanner.Scan() {
		f := strings.Split(scanner.Text(), "\t")
		if len(f) <= cols["AlternateAlleleVCF"] {
			continue
		}
		get := func(name string) string {
			idx := cols[name]
			if idx < 0 || idx >= len(f) {
				return ""
			}
			return f[idx]
		}

		if get("Assembly") != "GRCh38" {
			continue
		}
		variantType := classifyVariantType(get("Type"))
		if variantType == "" {
			continue
		}

		sig := get("ClinicalSignificance")
		if !isPathogenic(sig) {
			continue
		}

		name := get("Name")
		protein := extractProtein(name)
		if protein == "" {
			continue
		}

		posStr := get("PositionVCF")
		if posStr == "" || posStr == "-1" || strings.EqualFold(posStr, "na") {
			continue
		}
		pos, err := strconv.ParseInt(posStr, 10, 64)
		if err != nil || pos <= 0 {
			continue
		}

		ref := get("ReferenceAlleleVCF")
		alt := get("AlternateAlleleVCF")
		if ref == "" || alt == "" || strings.EqualFold(ref, "na") || strings.EqualFold(alt, "na") {
			continue
		}
		// Skip alleles that aren't real nucleotide sequences (ClinVar sometimes
		// stores "N/A" or symbolic values for structural variants).
		if !isNucleotideSeq(ref) || !isNucleotideSeq(alt) {
			continue
		}

		chrom := get("Chromosome")
		if chrom == "" || strings.EqualFold(chrom, "na") {
			continue
		}

		tx := extractTranscript(name)
		rev := ""
		if idx := cols["ReviewStatus"]; idx >= 0 && idx < len(f) {
			rev = f[idx]
		}

		entries = append(entries, SummaryEntry{
			Chrom:       chrom,
			Pos:         pos,
			Ref:         ref,
			Alt:         alt,
			Gene:        get("GeneSymbol"),
			Transcript:  tx,
			Protein:     protein,
			ClnSig:      sig,
			RevStatus:   rev,
			VariantType: variantType,
		})
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("clinvar summary scan: %w", err)
	}
	return entries, nil
}

// extractProtein returns the protein notation from a ClinVar Name field,
// e.g. "NM_000546.6(TP53):c.817C>T (p.Arg273Cys)" → "p.Arg273Cys".
// Returns "" if no protein annotation is present.
func extractProtein(name string) string {
	m := reProtein.FindStringSubmatch(name)
	if m == nil {
		return ""
	}
	p := "p." + m[1]
	// Skip entries where the protein change is just "p.?" or equals
	// p.Met1? (initiation codon uncertainty) or other non-specific values.
	if p == "p.?" || strings.HasSuffix(p, "?") {
		return ""
	}
	return p
}

// extractTranscript returns the unversioned RefSeq accession from a Name field,
// e.g. "NM_000546.6(TP53):c.817C>T (p.Arg273Cys)" → "NM_000546".
func extractTranscript(name string) string {
	m := reTranscript.FindStringSubmatch(name)
	if m == nil {
		return ""
	}
	// Strip version suffix (.6 etc.)
	tx := m[1]
	if dot := strings.IndexByte(tx, '.'); dot >= 0 {
		tx = tx[:dot]
	}
	return tx
}

// isPathogenic returns true for Pathogenic and Likely_pathogenic entries,
// but not "not Pathogenic" or conflicting interpretations.
func isPathogenic(sig string) bool {
	sig = strings.ToLower(sig)
	if strings.Contains(sig, "not pathogenic") {
		return false
	}
	if strings.Contains(sig, "conflicting") {
		return false
	}
	return strings.Contains(sig, "pathogenic")
}

// classifyVariantType maps ClinVar Type values to a normalized category.
// Returns "" for unsupported types (structural variants, microsatellites, etc.).
func classifyVariantType(clinvarType string) string {
	switch clinvarType {
	case "single nucleotide variant":
		return "snv"
	case "Deletion":
		return "deletion"
	case "Insertion":
		return "insertion"
	case "Indel":
		return "indel"
	}
	return ""
}

// isNucleotideSeq returns true if s consists only of A/C/G/T/N characters.
func isNucleotideSeq(s string) bool {
	if s == "" {
		return false
	}
	for _, c := range s {
		switch c {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			return false
		}
	}
	return true
}
func indexCols(header []string, names ...string) map[string]int {
	lookup := make(map[string]int, len(header))
	for i, h := range header {
		lookup[strings.TrimPrefix(h, "#")] = i
	}
	m := make(map[string]int, len(names))
	for _, n := range names {
		if idx, ok := lookup[n]; ok {
			m[n] = idx
		} else {
			m[n] = -1
		}
	}
	return m
}
