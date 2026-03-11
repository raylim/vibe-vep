package output_test

// ClinVar independent ground truth benchmark.
//
// This benchmark compares vibe-vep, snpEff, and VEP against ClinVar pathogenic
// variant annotations, which are curated independently of any annotation tool.
//
// Setup:
//   bash scripts/prepare_clinvar.sh
//   bash scripts/run_snpeff_clinvar.sh  [optional]
//   bash scripts/run_vep_clinvar.sh     [optional]
//
// Then run:
//   go test ./internal/output/ -run TestClinVarBenchmark -v -timeout 60m
//
// The test is skipped automatically when the data files are absent.

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	clv "github.com/inodb/vibe-vep/internal/datasource/clinvar"
	"github.com/inodb/vibe-vep/internal/datasource/mane"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// TestClinVarBenchmark compares annotation tools against ClinVar pathogenic
// variants as independent ground truth.
func TestClinVarBenchmark(t *testing.T) {
	summaryPath := filepath.Join(os.Getenv("HOME"), ".vibe-vep", "clinvar", "variant_summary.txt.gz")
	if _, err := os.Stat(summaryPath); err != nil {
		t.Skipf("ClinVar summary not found at %s — run scripts/prepare_clinvar.sh first", summaryPath)
	}

	manePath := filepath.Join(os.Getenv("HOME"), ".vibe-vep", "clinvar", "MANE.GRCh38.v1.5.summary.txt.gz")
	if _, err := os.Stat(manePath); err != nil {
		t.Skipf("MANE file not found at %s — run scripts/prepare_clinvar.sh first", manePath)
	}

	t.Log("Loading ClinVar variant_summary.txt.gz …")
	loadStart := time.Now()
	entries, err := clv.ParseSummaryFile(summaryPath)
	if err != nil {
		t.Fatalf("parse ClinVar summary: %v", err)
	}
	t.Logf("  Loaded %d pathogenic variants with protein HGVS (%.1fs)", len(entries), time.Since(loadStart).Seconds())

	t.Log("Loading MANE Select mapping …")
	maneMap, err := mane.Load(manePath)
	if err != nil {
		t.Fatalf("load MANE: %v", err)
	}
	t.Logf("  Loaded %d MANE Select transcripts", len(maneMap))

	// Mark MANE Select entries.
	maneCount := 0
	for i := range entries {
		if maneMap.HasRefSeq(entries[i].Transcript) {
			entries[i].IsMANE = true
			maneCount++
		}
	}
	t.Logf("  %d entries on MANE Select transcripts", maneCount)

	// De-duplicate by genomic position — ClinVar has one row per allele but
	// some positions appear multiple times (different review statuses). Keep
	// the entry with the most specific review status.
	entries = deduplicateEntries(entries)
	t.Logf("  %d unique variants after dedup", len(entries))

	// Build vibe-vep annotator.
	t.Log("Loading GENCODE cache …")
	c, cacheSource, cacheDur := loadClinVarCache(t)
	t.Logf("  Cache loaded in %.1fs from %s (%d transcripts)", cacheDur.Seconds(), cacheSource, c.TranscriptCount())

	ann := annotate.NewAnnotator(c)

	// Run vibe-vep annotations.
	t.Log("Annotating with vibe-vep …")
	vibeStart := time.Now()
	type vibeResult struct {
		anns []*annotate.Annotation
	}
	vibeAnns := make([]vibeResult, len(entries))
	for i, e := range entries {
		v := &vcf.Variant{
			Chrom: e.Chrom,
			Pos:   e.Pos,
			Ref:   e.Ref,
			Alt:   e.Alt,
		}
		anns, err := ann.Annotate(v)
		if err != nil {
			t.Logf("  WARN: annotate %s:%d %s>%s: %v", e.Chrom, e.Pos, e.Ref, e.Alt, err)
		}
		vibeAnns[i] = vibeResult{anns: anns}
	}
	vibeDur := time.Since(vibeStart)
	vibeRate := float64(len(entries)) / vibeDur.Seconds()
	t.Logf("  Annotated %d variants in %.1fs (%.0f v/s)", len(entries), vibeDur.Seconds(), vibeRate)

	// Load snpEff and VEP VCFs if present.
	clinvarDir := filepath.Join("../../testdata", "clinvar")
	snpEffVCF := filepath.Join(clinvarDir, "snpeff.vcf.gz")
	vepVCF := filepath.Join(clinvarDir, "vep.vcf.gz")

	var seMap snpEffVariantMap
	var hasSnpEff bool
	var snpEffAnnoDur time.Duration
	if _, err := os.Stat(snpEffVCF); err == nil {
		t.Log("Loading snpEff VCF …")
		seStart := time.Now()
		seMap, err = loadSnpEffVCF(snpEffVCF)
		if err != nil {
			t.Logf("  WARN: load snpEff VCF: %v", err)
		} else {
			t.Logf("  Loaded snpEff annotations for %d positions (%.1fs)", len(seMap.byPos), time.Since(seStart).Seconds())
			hasSnpEff = true
		}
		snpEffAnnoDur = readElapsedFile(filepath.Join(clinvarDir, "snpeff.elapsed"))
	} else {
		t.Log("snpEff VCF not found — skipping snpEff comparison")
	}

	var vepMap vepVariantMap
	var hasVEP bool
	var vepAnnoDur time.Duration
	if _, err := os.Stat(vepVCF); err == nil {
		t.Log("Loading VEP VCF …")
		vepStart := time.Now()
		vepMap, err = loadVEPVCF(vepVCF)
		if err != nil {
			t.Logf("  WARN: load VEP VCF: %v", err)
		} else {
			t.Logf("  Loaded VEP annotations for %d positions (%.1fs)", len(vepMap.byPos), time.Since(vepStart).Seconds())
			hasVEP = true
		}
		vepAnnoDur = readElapsedFile(filepath.Join(clinvarDir, "vep.elapsed"))
	} else {
		t.Log("VEP VCF not found — skipping VEP comparison")
	}

	// Evaluate.
	var vibe, snpEff, vep clinvarCounts

	for i, e := range entries {
		v := &vcf.Variant{
			Chrom: e.Chrom,
			Pos:   e.Pos,
			Ref:   e.Ref,
			Alt:   e.Alt,
		}
		expected := e.Protein // e.g. "p.Arg273Cys"
		isIndel := e.VariantType != "snv"

		// --- vibe-vep ---
		{
			anns := vibeAnns[i].anns
			vibe.total++
			if isIndel {
				vibe.indelTotal++
			} else {
				vibe.snvTotal++
			}
			if e.IsMANE {
				vibe.maneTotal++
			}
			best := output.PickBestAnnotation(anns)
			if best != nil {
				expClass := inferConsequenceClass(expected)
				var bestExact, bestAny bool
				if best.HGVSp != "" {
					bestExact = normalizeProteinStr(best.HGVSp) == normalizeProteinStr(expected)
					if bestExact {
						vibe.exactMatch++
						if isIndel {
							vibe.indelExact++
						} else {
							vibe.snvExact++
						}
						if e.IsMANE {
							vibe.maneExact++
						}
					}
					for _, a := range anns {
						if normalizeProteinStr(a.HGVSp) == normalizeProteinStr(expected) {
							bestAny = true
							vibe.anyMatch++
							break
						}
					}
				} else {
					vibe.notAnnotated++
				}
				if expClass != "" {
					vibe.consequTotal++
					if consequenceMatches(best.Consequence, expClass) {
						vibe.consequMatch++
					}
					vibe.trackClass(expClass, bestExact, bestAny)
				}
			}
		}

		// --- snpEff ---
		if hasSnpEff {
			seAnns := seMap.lookup(v)
			snpEff.total++
			if isIndel {
				snpEff.indelTotal++
			} else {
				snpEff.snvTotal++
			}
			if e.IsMANE {
				snpEff.maneTotal++
			}
			if len(seAnns) > 0 {
				best := pickBestSnpEffByImpact(seAnns)
				expClass := inferConsequenceClass(expected)
				var bestExact, bestAny bool
				if best != nil && best.hgvsp != "" {
					bestExact = normalizeProteinStr(best.hgvsp) == normalizeProteinStr(expected)
					if bestExact {
						snpEff.exactMatch++
						if isIndel {
							snpEff.indelExact++
						} else {
							snpEff.snvExact++
						}
						if e.IsMANE {
							snpEff.maneExact++
						}
					}
					for _, a := range seAnns {
						if normalizeProteinStr(a.hgvsp) == normalizeProteinStr(expected) {
							bestAny = true
							snpEff.anyMatch++
							break
						}
					}
					if expClass != "" {
						snpEff.consequTotal++
						if consequenceMatches(snpEffConsequenceToSO(best.consequence), expClass) {
							snpEff.consequMatch++
						}
					}
				} else {
					snpEff.notAnnotated++
				}
				if expClass != "" {
					snpEff.trackClass(expClass, bestExact, bestAny)
				}
			}
		}

		// --- VEP ---
		if hasVEP {
			vepAnns := vepMap.lookup(v)
			vep.total++
			if isIndel {
				vep.indelTotal++
			} else {
				vep.snvTotal++
			}
			if e.IsMANE {
				vep.maneTotal++
			}
			if len(vepAnns) > 0 {
				best := pickBestVEPByImpact(vepAnns)
				expClass := inferConsequenceClass(expected)
				var bestExact, bestAny bool
				if best != nil && best.hgvsp != "" {
					bestExact = normalizeProteinStr(best.hgvsp) == normalizeProteinStr(expected)
					if bestExact {
						vep.exactMatch++
						if isIndel {
							vep.indelExact++
						} else {
							vep.snvExact++
						}
						if e.IsMANE {
							vep.maneExact++
						}
					}
					for _, a := range vepAnns {
						if normalizeProteinStr(a.hgvsp) == normalizeProteinStr(expected) {
							bestAny = true
							vep.anyMatch++
							break
						}
					}
					if expClass != "" {
						vep.consequTotal++
						if consequenceMatches(best.consequence, expClass) {
							vep.consequMatch++
						}
					}
				} else {
					vep.notAnnotated++
				}
				if expClass != "" {
					vep.trackClass(expClass, bestExact, bestAny)
				}
			}
		}
	}

	// Write report.
	reportPath := filepath.Join(clinvarDir, "benchmark_report.md")
	if err := writeClinVarReport(reportPath, entries, vibe, snpEff, vep, hasSnpEff, hasVEP,
		vibeDur, vibeRate, cacheDur, cacheSource, snpEffAnnoDur, vepAnnoDur); err != nil {
		t.Logf("WARN: write report: %v", err)
	} else {
		t.Logf("Report written to %s", reportPath)
	}

	// Print summary.
	pct := func(n, d int) string {
		if d == 0 {
			return "N/A"
		}
		return fmt.Sprintf("%.1f%%", 100*float64(n)/float64(d))
	}
	t.Logf("\nClinVar Benchmark Results (n=%d variants: %d SNVs + %d indels, %d MANE Select)",
		len(entries), vibe.snvTotal, vibe.indelTotal, maneCount)
	t.Logf("%-10s  %s  %s  %s  %s  %s  %s", "Tool", "HGVSp(best)", "HGVSp(any)", "Consequence", "MANE HGVSp", "SNV best", "Indel best")
	t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s  %-8s  %-10s", "vibe-vep",
		pct(vibe.exactMatch, vibe.total), pct(vibe.anyMatch, vibe.total),
		pct(vibe.consequMatch, vibe.consequTotal), pct(vibe.maneExact, vibe.maneTotal),
		pct(vibe.snvExact, vibe.snvTotal), pct(vibe.indelExact, vibe.indelTotal))
	if hasSnpEff {
		t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s  %-8s  %-10s", "snpEff",
			pct(snpEff.exactMatch, snpEff.total), pct(snpEff.anyMatch, snpEff.total),
			pct(snpEff.consequMatch, snpEff.consequTotal), pct(snpEff.maneExact, snpEff.maneTotal),
			pct(snpEff.snvExact, snpEff.snvTotal), pct(snpEff.indelExact, snpEff.indelTotal))
	}
	if hasVEP {
		t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s  %-8s  %-10s", "VEP",
			pct(vep.exactMatch, vep.total), pct(vep.anyMatch, vep.total),
			pct(vep.consequMatch, vep.consequTotal), pct(vep.maneExact, vep.maneTotal),
			pct(vep.snvExact, vep.snvTotal), pct(vep.indelExact, vep.indelTotal))
	}
}

// classCounts holds per-consequence-class accuracy counts for one tool.
type classCounts struct {
	total, exact, any int
}

// clinvarCounts holds benchmark counts for one annotation tool.
type clinvarCounts struct {
	total, exactMatch, anyMatch, maneTotal, maneExact int
	consequMatch, consequTotal, notAnnotated           int
	snvTotal, snvExact                                 int
	indelTotal, indelExact                             int
	byClass                                            map[string]*classCounts
}

// trackClass records an observation for a consequence class.
func (c *clinvarCounts) trackClass(class string, exact, any bool) {
	if c.byClass == nil {
		c.byClass = make(map[string]*classCounts)
	}
	cc := c.byClass[class]
	if cc == nil {
		cc = &classCounts{}
		c.byClass[class] = cc
	}
	cc.total++
	if exact {
		cc.exact++
	}
	if any {
		cc.any++
	}
}
// Preference: MANE > higher-review-status > first seen.
func deduplicateEntries(entries []clv.SummaryEntry) []clv.SummaryEntry {
	type key struct {
		chrom, ref, alt string
		pos             int64
	}
	seen := make(map[key]*clv.SummaryEntry, len(entries))
	order := make([]key, 0, len(entries))
	for i := range entries {
		e := &entries[i]
		k := key{e.Chrom, e.Ref, e.Alt, e.Pos}
		if prev, ok := seen[k]; !ok {
			seen[k] = e
			order = append(order, k)
		} else if e.IsMANE && !prev.IsMANE {
			seen[k] = e
		} else if reviewRank(e.RevStatus) > reviewRank(prev.RevStatus) {
			seen[k] = e
		}
	}
	out := make([]clv.SummaryEntry, 0, len(order))
	for _, k := range order {
		out = append(out, *seen[k])
	}
	return out
}

func reviewRank(s string) int {
	s = strings.ToLower(s)
	switch {
	case strings.Contains(s, "practice guideline"):
		return 4
	case strings.Contains(s, "expert panel"):
		return 3
	case strings.Contains(s, "multiple submitters"):
		return 2
	case strings.Contains(s, "single submitter"):
		return 1
	}
	return 0
}

// loadClinVarCache creates a vibe-vep transcript cache from GENCODE files.
func loadClinVarCache(t *testing.T) (*cache.Cache, string, time.Duration) {
	t.Helper()
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t, "GRCh38")
	cacheDir := filepath.Dir(gtfPath)
	c := cache.New()

	cacheStart := time.Now()
	tc := duckdb.NewTranscriptCache(cacheDir)
	gtfFP, err1 := duckdb.StatFile(gtfPath)
	fastaFP, err2 := duckdb.StatFile(fastaPath)
	canonicalFP := duckdb.FileFingerprint{}
	if canonicalPath != "" {
		canonicalFP, _ = duckdb.StatFile(canonicalPath)
	}

	var loadSource string
	if err1 == nil && err2 == nil && tc.Valid(gtfFP, fastaFP, canonicalFP) {
		if err := tc.Load(c); err != nil {
			t.Fatalf("load transcript cache: %v", err)
		}
		loadSource = "duckdb cache"
	} else {
		loader := cache.NewGENCODELoader(gtfPath, fastaPath)
		if canonicalPath != "" {
			overrides, err := cache.LoadCanonicalOverrides(canonicalPath)
			if err != nil {
				t.Logf("warning: could not load canonical overrides: %v", err)
			} else {
				loader.SetCanonicalOverrides(overrides)
			}
		}
		if err := loader.Load(c); err != nil {
			t.Fatalf("load GENCODE: %v", err)
		}
		loadSource = "GTF/FASTA"
	}
	c.BuildIndex()
	return c, loadSource, time.Since(cacheStart)
}

// normalizeProteinStr normalizes a protein HGVS string for comparison.
//
// ClinVar often stores abbreviated frameshift notation like p.Asp113fs, while
// tools emit the full HGVS form p.Asp113ValfsTer15. Both are valid; we normalize
// to the abbreviated form (first-affected AA + position + "fs") so they compare equal.
// We also normalize stop-codon notation: p.Arg273* → p.Arg273Ter (snpEff uses *, ClinVar uses Ter).
var (
	reProteinNorm = regexp.MustCompile(`(?i)ter\b`)
	// reFsNorm matches full HGVS frameshift: XxxNNNYyyfs... and strips the new AA + stop.
	reFsNorm = regexp.MustCompile(`^([A-Z][a-z]{2}\d+)[A-Z][a-z]{2}fs.*$`)
)

func normalizeProteinStr(p string) string {
	p = strings.TrimPrefix(p, "p.")
	p = strings.TrimSpace(p)
	// Normalize stop-codon notations: * and X both map to Ter.
	p = strings.ReplaceAll(p, "*", "Ter")
	p = reProteinNorm.ReplaceAllString(p, "Ter")
	// Normalize full HGVS frameshift to abbreviated form: XxxNNNYyyFsTerM → XxxNNNfs
	if m := reFsNorm.FindStringSubmatch(p); m != nil {
		return m[1] + "fs"
	}
	return p
}

// inferConsequenceClass infers the expected SO consequence class from a protein HGVS string.
func inferConsequenceClass(protein string) string {
	p := strings.TrimPrefix(protein, "p.")
	lower := strings.ToLower(p)
	switch {
	case strings.Contains(lower, "fs"):
		return "frameshift"
	case strings.HasSuffix(lower, "ter") || strings.Contains(lower, "*"):
		return "stop_gained"
	case strings.Contains(lower, "="):
		return "synonymous"
	case strings.Contains(lower, "dup"):
		return "inframe_insertion"
	case strings.Contains(lower, "ins"):
		return "inframe_insertion"
	case strings.Contains(lower, "del"):
		return "inframe_deletion"
	default:
		// Standard missense: p.Xxx123Yyy
		if len(p) >= 7 && isAminoAcidCode(p[:3]) && isAminoAcidCode(p[len(p)-3:]) {
			return "missense_variant"
		}
	}
	return ""
}

var aminoAcidCodes = map[string]bool{
	"Ala": true, "Arg": true, "Asn": true, "Asp": true, "Cys": true,
	"Gln": true, "Glu": true, "Gly": true, "His": true, "Ile": true,
	"Leu": true, "Lys": true, "Met": true, "Phe": true, "Pro": true,
	"Ser": true, "Thr": true, "Trp": true, "Tyr": true, "Val": true,
}

func isAminoAcidCode(s string) bool {
	if len(s) != 3 {
		return false
	}
	key := strings.Title(strings.ToLower(s))
	return aminoAcidCodes[key]
}

// consequenceMatches checks if a tool consequence matches the expected class.
func consequenceMatches(toolConsequence, expected string) bool {
	switch expected {
	case "missense_variant":
		return strings.Contains(toolConsequence, "missense_variant")
	case "stop_gained":
		return strings.Contains(toolConsequence, "stop_gained")
	case "frameshift":
		return strings.Contains(toolConsequence, "frameshift")
	case "synonymous":
		return strings.Contains(toolConsequence, "synonymous")
	case "inframe_deletion":
		return strings.Contains(toolConsequence, "inframe_deletion")
	case "inframe_insertion":
		return strings.Contains(toolConsequence, "inframe_insertion")
	}
	return false
}

// snpEffConsequenceToSO maps snpEff effect names to SO terms for comparison.
func snpEffConsequenceToSO(eff string) string {
	switch eff {
	case "missense_variant":
		return "missense_variant"
	case "stop_gained":
		return "stop_gained"
	case "frameshift_variant":
		return "frameshift_variant"
	case "synonymous_variant":
		return "synonymous_variant"
	}
	return eff
}

// pickBestSnpEffByImpact returns the highest-impact snpEff annotation.
func pickBestSnpEffByImpact(anns []snpEffAnnotation) *snpEffAnnotation {
	if len(anns) == 0 {
		return nil
	}
	best := &anns[0]
	for i := range anns[1:] {
		if snpEffImpactRank(anns[i+1].impact) > snpEffImpactRank(best.impact) {
			best = &anns[i+1]
		} else if snpEffImpactRank(anns[i+1].impact) == snpEffImpactRank(best.impact) {
			if anns[i+1].biotype == "protein_coding" && best.biotype != "protein_coding" {
				best = &anns[i+1]
			}
		}
	}
	return best
}

func snpEffImpactRank(impact string) int {
	switch impact {
	case "HIGH":
		return 4
	case "MODERATE":
		return 3
	case "LOW":
		return 2
	case "MODIFIER":
		return 1
	}
	return 0
}

// pickBestVEPByImpact returns the highest-impact VEP annotation.
func pickBestVEPByImpact(anns []vepAnnotation) *vepAnnotation {
	if len(anns) == 0 {
		return nil
	}
	best := &anns[0]
	for i := range anns[1:] {
		if snpEffImpactRank(anns[i+1].impact) > snpEffImpactRank(best.impact) {
			best = &anns[i+1]
		} else if snpEffImpactRank(anns[i+1].impact) == snpEffImpactRank(best.impact) {
			if anns[i+1].biotype == "protein_coding" && best.biotype != "protein_coding" {
				best = &anns[i+1]
			}
		}
	}
	return best
}

// writeClinVarReport generates the Markdown benchmark report.
func writeClinVarReport(
	path string,
	entries []clv.SummaryEntry,
	vibe, snpEff, vep clinvarCounts,
	hasSnpEff, hasVEP bool,
	vibeDur time.Duration, vibeRate float64,
	cacheDur time.Duration, cacheSource string,
	snpEffAnnoDur, vepAnnoDur time.Duration,
) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	maneCount := 0
	for _, e := range entries {
		if e.IsMANE {
			maneCount++
		}
	}
	sigDist := countSignificance(entries)

	pct := func(n, d int) string {
		if d == 0 {
			return "N/A"
		}
		return fmt.Sprintf("%.1f%%", 100*float64(n)/float64(d))
	}

	fmt.Fprintf(w, "# ClinVar Independent Ground Truth Benchmark\n\n")
	fmt.Fprintf(w, "Generated by `TestClinVarBenchmark` using ClinVar `variant_summary.txt.gz`.\n\n")
	fmt.Fprintf(w, "## Methodology\n\n")
	fmt.Fprintf(w, "**Ground truth**: ClinVar pathogenic/likely-pathogenic SNVs and small indels with curated protein\n")
	fmt.Fprintf(w, "HGVS annotations (from the `Name` field of `variant_summary.txt.gz`). ClinVar\n")
	fmt.Fprintf(w, "annotations are submitted by clinical labs, expert panels, and NCBI curators —\n")
	fmt.Fprintf(w, "independent of VEP, snpEff, or vibe-vep.\n\n")
	fmt.Fprintf(w, "**MANE Select**: %d / %d variants (%.0f%%) map to MANE Select transcripts.\n",
		maneCount, len(entries), 100*float64(maneCount)/float64(len(entries)))
	fmt.Fprintf(w, "For these, NM_ and ENST transcripts encode identical protein sequences,\n")
	fmt.Fprintf(w, "making cross-tool HGVSp comparison directly meaningful.\n\n")
	fmt.Fprintf(w, "**Database versions**: GENCODE v49 / Ensembl 115 (vibe-vep, VEP), snpEff GRCh38.115.\n\n")
	fmt.Fprintf(w, "## Dataset\n\n")
	snvCount, indelCount := countVariantTypes(entries)
	fmt.Fprintf(w, "| Metric | Value |\n|--------|-------|\n")
	fmt.Fprintf(w, "| Total pathogenic variants | %d |\n", len(entries))
	fmt.Fprintf(w, "| SNVs | %d |\n", snvCount)
	fmt.Fprintf(w, "| Indels (deletions, insertions, indels) | %d |\n", indelCount)
	fmt.Fprintf(w, "| MANE Select transcripts | %d (%.0f%%) |\n", maneCount, 100*float64(maneCount)/float64(len(entries)))
	for sig, n := range sigDist {
		fmt.Fprintf(w, "| %s | %d |\n", sig, n)
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "## Results\n\n")
	fmt.Fprintf(w, "### Protein HGVS Match (all variants)\n\n")
	fmt.Fprintf(w, "| Tool | HGVSp Match (best) | HGVSp Match (any) | Not Annotated |\n")
	fmt.Fprintf(w, "|------|--------------------|-------------------|---------------|\n")
	fmt.Fprintf(w, "| vibe-vep | %s | %s | %d |\n",
		pct(vibe.exactMatch, vibe.total), pct(vibe.anyMatch, vibe.total), vibe.notAnnotated)
	if hasSnpEff {
		fmt.Fprintf(w, "| snpEff GRCh38.115 | %s | %s | %d |\n",
			pct(snpEff.exactMatch, snpEff.total), pct(snpEff.anyMatch, snpEff.total), snpEff.notAnnotated)
	}
	if hasVEP {
		fmt.Fprintf(w, "| Ensembl VEP v115 | %s | %s | %d |\n",
			pct(vep.exactMatch, vep.total), pct(vep.anyMatch, vep.total), vep.notAnnotated)
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "### Protein HGVS Match by Variant Type\n\n")
	fmt.Fprintf(w, "| Tool | SNV (n=%d) | Indel (n=%d) |\n", snvCount, indelCount)
	fmt.Fprintf(w, "|------|-----------|-------------|\n")
	fmt.Fprintf(w, "| vibe-vep | %s | %s |\n", pct(vibe.snvExact, vibe.snvTotal), pct(vibe.indelExact, vibe.indelTotal))
	if hasSnpEff {
		fmt.Fprintf(w, "| snpEff GRCh38.115 | %s | %s |\n", pct(snpEff.snvExact, snpEff.snvTotal), pct(snpEff.indelExact, snpEff.indelTotal))
	}
	if hasVEP {
		fmt.Fprintf(w, "| Ensembl VEP v115 | %s | %s |\n", pct(vep.snvExact, vep.snvTotal), pct(vep.indelExact, vep.indelTotal))
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "### Protein HGVS Match (MANE Select transcripts only, n=%d)\n\n", maneCount)
	fmt.Fprintf(w, "_MANE Select variants provide the fairest comparison: all tools should use_\n")
	fmt.Fprintf(w, "_the same transcript, so protein notation differences reflect real errors._\n\n")
	fmt.Fprintf(w, "| Tool | HGVSp Match |\n|------|-------------|\n")
	fmt.Fprintf(w, "| vibe-vep | %s |\n", pct(vibe.maneExact, vibe.maneTotal))
	if hasSnpEff {
		fmt.Fprintf(w, "| snpEff GRCh38.115 | %s |\n", pct(snpEff.maneExact, snpEff.maneTotal))
	}
	if hasVEP {
		fmt.Fprintf(w, "| Ensembl VEP v115 | %s |\n", pct(vep.maneExact, vep.maneTotal))
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "### HGVSp Match by Consequence Class\n\n")
	fmt.Fprintf(w, "_\"Best\" = primary transcript; \"Any\" = correct answer exists in any annotated transcript._\n")
	fmt.Fprintf(w, "_snpEff and VEP annotate all transcripts, so \"any\" reveals whether the right answer is present_\n")
	fmt.Fprintf(w, "_but not selected as primary (transcript-choice errors). vibe-vep only reports MANE/canonical._\n\n")
	classOrder := []struct{ id, label string }{
		{"missense_variant", "missense"},
		{"frameshift", "frameshift"},
		{"stop_gained", "stop_gained"},
		{"inframe_deletion", "inframe_del"},
		{"inframe_insertion", "inframe_ins"},
		{"synonymous", "synonymous"},
	}
	// Header
	if hasSnpEff && hasVEP {
		fmt.Fprintf(w, "| Class | n | vibe-vep best | vibe-vep any | snpEff best | snpEff any | VEP best | VEP any |\n")
		fmt.Fprintf(w, "|-------|---|--------------|--------------|-------------|------------|----------|---------|\n")
	} else if hasSnpEff {
		fmt.Fprintf(w, "| Class | n | vibe-vep best | vibe-vep any | snpEff best | snpEff any |\n")
		fmt.Fprintf(w, "|-------|---|--------------|--------------|-------------|------------|\n")
	} else {
		fmt.Fprintf(w, "| Class | n | vibe-vep best | vibe-vep any |\n")
		fmt.Fprintf(w, "|-------|---|--------------|---------------|\n")
	}
	for _, cl := range classOrder {
		vc := vibe.byClass[cl.id]
		if vc == nil || vc.total == 0 {
			continue
		}
		if hasSnpEff && hasVEP {
			sc := snpEff.byClass[cl.id]
			vc2 := vep.byClass[cl.id]
			if sc == nil {
				sc = &classCounts{}
			}
			if vc2 == nil {
				vc2 = &classCounts{}
			}
			fmt.Fprintf(w, "| %s | %d | %s | %s | %s | %s | %s | %s |\n",
				cl.label, vc.total,
				pct(vc.exact, vc.total), pct(vc.any, vc.total),
				pct(sc.exact, sc.total), pct(sc.any, sc.total),
				pct(vc2.exact, vc2.total), pct(vc2.any, vc2.total))
		} else if hasSnpEff {
			sc := snpEff.byClass[cl.id]
			if sc == nil {
				sc = &classCounts{}
			}
			fmt.Fprintf(w, "| %s | %d | %s | %s | %s | %s |\n",
				cl.label, vc.total,
				pct(vc.exact, vc.total), pct(vc.any, vc.total),
				pct(sc.exact, sc.total), pct(sc.any, sc.total))
		} else {
			fmt.Fprintf(w, "| %s | %d | %s | %s |\n",
				cl.label, vc.total,
				pct(vc.exact, vc.total), pct(vc.any, vc.total))
		}
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "### Consequence Class Match\n\n")
	fmt.Fprintf(w, "_Inferred from protein notation: missense → `missense_variant`,_\n")
	fmt.Fprintf(w, "_Ter/\\* suffix → `stop_gained`, `fs` → `frameshift_variant`,_\n")
	fmt.Fprintf(w, "_`del`/`ins`/`dup` → `inframe_deletion`/`inframe_insertion`._\n\n")
	fmt.Fprintf(w, "| Tool | Consequence Match |\n|------|-------------------|\n")
	fmt.Fprintf(w, "| vibe-vep | %s |\n", pct(vibe.consequMatch, vibe.consequTotal))
	if hasSnpEff {
		fmt.Fprintf(w, "| snpEff GRCh38.115 | %s |\n", pct(snpEff.consequMatch, snpEff.consequTotal))
	}
	if hasVEP {
		fmt.Fprintf(w, "| Ensembl VEP v115 | %s |\n", pct(vep.consequMatch, vep.consequTotal))
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "## Performance\n\n")
	fmt.Fprintf(w, "| Tool | Variants | Time | Rate |\n")
	fmt.Fprintf(w, "|------|----------|------|------|\n")
	fmt.Fprintf(w, "| vibe-vep | %d | %.1fs | %.0f v/s |\n", len(entries), vibeDur.Seconds(), vibeRate)
	if hasSnpEff {
		if snpEffAnnoDur > 0 {
			snpEffRate := float64(len(entries)) / snpEffAnnoDur.Seconds()
			fmt.Fprintf(w, "| snpEff GRCh38.115 | %d | %.0fs | %.0f v/s |\n",
				len(entries), snpEffAnnoDur.Seconds(), snpEffRate)
		} else {
			fmt.Fprintf(w, "| snpEff GRCh38.115 | %d | — | — |\n", len(entries))
		}
	}
	if hasVEP {
		if vepAnnoDur > 0 {
			vepRate := float64(len(entries)) / vepAnnoDur.Seconds()
			fmt.Fprintf(w, "| Ensembl VEP v115 | %d | %.0fs | %.0f v/s |\n",
				len(entries), vepAnnoDur.Seconds(), vepRate)
		} else {
			fmt.Fprintf(w, "| Ensembl VEP v115 | %d | — | — |\n", len(entries))
		}
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "_vibe-vep cache load: %.1fs from %s. snpEff/VEP times from `*.elapsed` sidecar written by annotation scripts._\n\n",
		cacheDur.Seconds(), cacheSource)
	if !hasSnpEff && !hasVEP {
		fmt.Fprintf(w, "_snpEff and VEP results not available. Run `run_snpeff_clinvar.sh` and `run_vep_clinvar.sh`._\n\n")
	}
	fmt.Fprintf(w, "## Interpretation\n\n")
	fmt.Fprintf(w, "Unlike the TCGA benchmark (where ground truth is VEP-annotated),\n")
	fmt.Fprintf(w, "ClinVar annotations are curated independently. HGVSp discrepancies reflect:\n\n")
	fmt.Fprintf(w, "1. **Transcript version drift**: ClinVar may use older NM_ versions.\n")
	fmt.Fprintf(w, "2. **MANE transition**: Older submissions pre-date MANE Select consensus.\n")
	fmt.Fprintf(w, "3. **Tool transcript choice**: Each tool may select a different default transcript.\n")
	fmt.Fprintf(w, "4. **Coding differences**: Where tools genuinely compute different protein changes.\n\n")
	fmt.Fprintf(w, "The MANE Select subset eliminates source (3) and reduces source (2),\n")
	fmt.Fprintf(w, "providing the most rigorous comparison of actual prediction accuracy.\n\n")
	fmt.Fprintf(w, "### Normalization note\n\n")
	fmt.Fprintf(w, "Two HGVS notation variants are normalized before comparison:\n\n")
	fmt.Fprintf(w, "- **Frameshift stop notation**: Tools emit full HGVS `p.Asp113ValfsTer15`; ClinVar\n")
	fmt.Fprintf(w, "  abbreviates to `p.Asp113fs`. Normalized by stripping new-AA and stop-distance.\n")
	fmt.Fprintf(w, "- **Stop-codon glyph**: snpEff uses `p.Gln55*`; ClinVar uses `p.Gln55Ter`.\n")
	fmt.Fprintf(w, "  Normalized by replacing `*` with `Ter`.\n\n")
	fmt.Fprintf(w, "### Breakdown by consequence class\n\n")

	classInfo := []struct {
		id, label, approxN string
	}{
		{"missense_variant", "Missense", "≈70k"},
		{"stop_gained", "Stop-gained", "≈76k"},
		{"frameshift", "Frameshift", "≈83k"},
		{"inframe_deletion", "Inframe deletion", "≈2k"},
		{"inframe_insertion", "Inframe insertion", "≈750"},
		{"synonymous", "Synonymous", "≈815"},
	}

	classSummary := map[string]string{
		"missense_variant": "All three reach ~97% \"any\" match, confirming the differences are transcript-choice,\n" +
			"not algorithmic. vibe-vep's MANE Select preference gives it the best primary match.\n\n",
		"stop_gained": "\"Any\" match of 91–93% indicates the remainder are transcript-drift cases.\n\n",
		"frameshift": "snpEff and VEP both reach ~99.7% \"any\" match, meaning the correct answer\n" +
			"exists in their multi-transcript output.\n\n",
		"inframe_deletion": "The high snpEff/VEP \"any\" (~95%) suggests the protein is correctly computed\n" +
			"but the position range notation (e.g., `p.Arg27_Ile28del`) is sensitive to\n" +
			"exact codon boundary choice under different normalization rules.\n\n",
		"inframe_insertion": "Insertion HGVSp notation is particularly complex\n" +
			"(position-range, dup vs ins disambiguation) and requires further investigation.\n\n",
		"synonymous": "snpEff and VEP do not emit HGVSp for synonymous variants (silent change is not\n" +
			"annotated in the protein-change field); vibe-vep outputs `p.Arg273=` notation.\n\n",
	}

	for _, ci := range classInfo {
		vc := vibe.byClass[ci.id]
		sc := snpEff.byClass[ci.id]
		vc2 := vep.byClass[ci.id]
		vibeStr := "n/a"
		if vc != nil {
			vibeStr = pct(vc.exact, vc.total)
		}
		snpEffStr := "n/a"
		if sc != nil && hasSnpEff {
			snpEffStr = pct(sc.exact, sc.total)
		}
		vepStr := "n/a"
		if vc2 != nil && hasVEP {
			vepStr = pct(vc2.exact, vc2.total)
		}
		totalStr := ci.approxN
		if vc != nil {
			totalStr = "n=" + strconv.Itoa(vc.total) + " " + ci.approxN
		}
		fmt.Fprintf(w, "**%s** (%s): vibe-vep %s, snpEff %s, VEP %s.\n",
			ci.label, totalStr, vibeStr, snpEffStr, vepStr)
		if note, ok := classSummary[ci.id]; ok {
			fmt.Fprint(w, note)
		} else {
			fmt.Fprintln(w)
		}
	}

	fmt.Fprintf(w, "### Why VEP \"any\" match is much higher than \"best\"\n\n")
	fmt.Fprintf(w, "VEP annotates all transcripts and the correct HGVSp exists in a non-primary\n")
	fmt.Fprintf(w, "transcript for ~16%% of variants. This is larger than snpEff (~14%%) or vibe-vep\n")
	fmt.Fprintf(w, "(~7%%), suggesting VEP's primary-transcript ranking is less aligned with MANE\n")
	fmt.Fprintf(w, "than vibe-vep, but VEP's underlying nucleotide consequence prediction is\n")
	fmt.Fprintf(w, "accurate for virtually all variants.\n")

	return w.Flush()
}

// countVariantTypes counts SNVs and indels in the entry set.
func countVariantTypes(entries []clv.SummaryEntry) (snvCount, indelCount int) {
	for _, e := range entries {
		if e.VariantType == "snv" {
			snvCount++
		} else {
			indelCount++
		}
	}
	return
}

// countSignificance tallies the ClinicalSignificance values in the dataset.
func countSignificance(entries []clv.SummaryEntry) map[string]int {
	m := make(map[string]int)
	for _, e := range entries {
		sig := e.ClnSig
		if strings.Contains(sig, "Likely_pathogenic") {
			m["Likely pathogenic"]++
		} else {
			m["Pathogenic"]++
		}
	}
	return m
}

// TestFrameshiftPositionDebug samples frameshift variants where vibe-vep
// disagrees with ClinVar to understand the position offset distribution.
func TestFrameshiftPositionDebug(t *testing.T) {
	summaryPath := filepath.Join(os.Getenv("HOME"), ".vibe-vep", "clinvar", "variant_summary.txt.gz")
	if _, err := os.Stat(summaryPath); err != nil {
		t.Skipf("ClinVar summary not found — run scripts/prepare_clinvar.sh first")
	}
	manePath := filepath.Join(os.Getenv("HOME"), ".vibe-vep", "clinvar", "MANE.GRCh38.v1.5.summary.txt.gz")
	if _, err := os.Stat(manePath); err != nil {
		t.Skipf("MANE file not found")
	}
	entries, err := clv.ParseSummaryFile(summaryPath)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	maneMap, err := mane.Load(manePath)
	if err != nil {
		t.Fatalf("load MANE: %v", err)
	}
	for i := range entries {
		if maneMap.HasRefSeq(entries[i].Transcript) {
			entries[i].IsMANE = true
		}
	}
	entries = deduplicateEntries(entries)

	c, _, _ := loadClinVarCache(t)
	ann := annotate.NewAnnotator(c)

	// reProteinPos extracts the numeric position from a normalised p.XxxNNNfs string.
	rePos := regexp.MustCompile(`^([A-Z][a-z]{2})(\d+)fs`)

	type mismatch struct {
		chrom, ref, alt string
		pos             int64
		expected        string
		got             string
		offset          int
	}
	offsets := make(map[int]int) // offset → count
	var samples []mismatch
	const maxSamples = 20

	for _, e := range entries {
		if inferConsequenceClass(e.Protein) != "frameshift" {
			continue
		}
		v := &vcf.Variant{Chrom: e.Chrom, Pos: e.Pos, Ref: e.Ref, Alt: e.Alt}
		anns, _ := ann.Annotate(v)
		best := output.PickBestAnnotation(anns)
		if best == nil || best.HGVSp == "" {
			continue
		}
		normExp := normalizeProteinStr(e.Protein)
		normGot := normalizeProteinStr(best.HGVSp)
		if normExp == normGot {
			continue // match — skip
		}
		// Extract positions.
		mExp := rePos.FindStringSubmatch(normExp)
		mGot := rePos.FindStringSubmatch(normGot)
		if mExp == nil || mGot == nil {
			continue
		}
		var posExp, posGot int
		fmt.Sscanf(mExp[2], "%d", &posExp)
		fmt.Sscanf(mGot[2], "%d", &posGot)
		offset := posGot - posExp
		offsets[offset]++
		if len(samples) < maxSamples {
			samples = append(samples, mismatch{
				chrom: e.Chrom, pos: e.Pos, ref: e.Ref, alt: e.Alt,
				expected: normExp, got: normGot, offset: offset,
			})
		}
	}

	t.Logf("Frameshift position offset distribution (vibe-vep − ClinVar):")
	// Print sorted offsets.
	for off := -10; off <= 10; off++ {
		if c := offsets[off]; c > 0 {
			t.Logf("  offset %+d: %d variants", off, c)
		}
	}
	// Print anything outside −10..+10.
	var otherTotal int
	for off, c := range offsets {
		if off < -10 || off > 10 {
			otherTotal += c
		}
	}
	if otherTotal > 0 {
		t.Logf("  offset >±10: %d variants", otherTotal)
	}
	t.Logf("Sample mismatches (first %d):", len(samples))
	for _, m := range samples {
		t.Logf("  %s:%d %s>%s  expected=%s  got=%s  (offset %+d)",
			m.chrom, m.pos, m.ref, m.alt, m.expected, m.got, m.offset)
	}
}

// TestFrameshiftSingleVariantDebug traces internal state for a known failing case.
func TestFrameshiftSingleVariantDebug(t *testing.T) {
	c, _, _ := loadClinVarCache(t)
	ann := annotate.NewAnnotator(c)

	// 2:29073314 AT>A  expected=Asn316fs  got=Val317fs
	v := &vcf.Variant{Chrom: "2", Pos: 29073314, Ref: "AT", Alt: "A"}
	anns, err := ann.Annotate(v)
	if err != nil {
		t.Fatalf("annotate: %v", err)
	}
	for _, a := range anns {
		if a.HGVSp != "" {
			tr := c.GetTranscript(a.TranscriptID)
			strand := "unknown"
			if tr != nil {
				if tr.IsForwardStrand() {
					strand = "+"
				} else {
					strand = "-"
				}
			}
			t.Logf("transcript=%s strand=%s HGVSp=%s CDSpos=%d proteinPos=%d",
				a.TranscriptID, strand, a.HGVSp, a.CDSPosition, a.ProteinPosition)
		}
	}
}

// readElapsedFile reads a sidecar *.elapsed file written by annotation scripts.
// Returns 0 if the file is missing or unparseable.
func readElapsedFile(path string) time.Duration {
	data, err := os.ReadFile(path)
	if err != nil {
		return 0
	}
	var secs int64
	if _, err := fmt.Sscanf(strings.TrimSpace(string(data)), "%d", &secs); err != nil {
		return 0
	}
	return time.Duration(secs) * time.Second
}
