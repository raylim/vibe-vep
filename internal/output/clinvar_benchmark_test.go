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
	t.Logf("  Loaded %d pathogenic SNVs with protein HGVS (%.1fs)", len(entries), time.Since(loadStart).Seconds())

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
	if _, err := os.Stat(snpEffVCF); err == nil {
		t.Log("Loading snpEff VCF …")
		seMap, err = loadSnpEffVCF(snpEffVCF)
		if err != nil {
			t.Logf("  WARN: load snpEff VCF: %v", err)
		} else {
			t.Logf("  Loaded snpEff annotations for %d positions", len(seMap.byPos))
			hasSnpEff = true
		}
	} else {
		t.Log("snpEff VCF not found — skipping snpEff comparison")
	}

	var vepMap vepVariantMap
	var hasVEP bool
	if _, err := os.Stat(vepVCF); err == nil {
		t.Log("Loading VEP VCF …")
		vepMap, err = loadVEPVCF(vepVCF)
		if err != nil {
			t.Logf("  WARN: load VEP VCF: %v", err)
		} else {
			t.Logf("  Loaded VEP annotations for %d positions", len(vepMap.byPos))
			hasVEP = true
		}
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

		// --- vibe-vep ---
		{
			anns := vibeAnns[i].anns
			vibe.total++
			if e.IsMANE {
				vibe.maneTotal++
			}
			best := output.PickBestAnnotation(anns)
			if best != nil {
				if best.HGVSp != "" {
					if normalizeProteinStr(best.HGVSp) == normalizeProteinStr(expected) {
						vibe.exactMatch++
						if e.IsMANE {
							vibe.maneExact++
						}
					}
					// Check any annotation.
					for _, a := range anns {
						if normalizeProteinStr(a.HGVSp) == normalizeProteinStr(expected) {
							vibe.anyMatch++
							break
						}
					}
				} else {
					vibe.notAnnotated++
				}
				// Consequence match: infer expected class from protein.
				expClass := inferConsequenceClass(expected)
				if expClass != "" {
					vibe.consequTotal++
					if consequenceMatches(best.Consequence, expClass) {
						vibe.consequMatch++
					}
				}
			}
		}

		// --- snpEff ---
		if hasSnpEff {
			seAnns := seMap.lookup(v)
			snpEff.total++
			if e.IsMANE {
				snpEff.maneTotal++
			}
			if len(seAnns) > 0 {
				best := pickBestSnpEffByImpact(seAnns)
				if best != nil && best.hgvsp != "" {
					if normalizeProteinStr(best.hgvsp) == normalizeProteinStr(expected) {
						snpEff.exactMatch++
						if e.IsMANE {
							snpEff.maneExact++
						}
					}
					for _, a := range seAnns {
						if normalizeProteinStr(a.hgvsp) == normalizeProteinStr(expected) {
							snpEff.anyMatch++
							break
						}
					}
					expClass := inferConsequenceClass(expected)
					if expClass != "" {
						snpEff.consequTotal++
						if consequenceMatches(snpEffConsequenceToSO(best.consequence), expClass) {
							snpEff.consequMatch++
						}
					}
				} else {
					snpEff.notAnnotated++
				}
			}
		}

		// --- VEP ---
		if hasVEP {
			vepAnns := vepMap.lookup(v)
			vep.total++
			if e.IsMANE {
				vep.maneTotal++
			}
			if len(vepAnns) > 0 {
				best := pickBestVEPByImpact(vepAnns)
				if best != nil && best.hgvsp != "" {
					if normalizeProteinStr(best.hgvsp) == normalizeProteinStr(expected) {
						vep.exactMatch++
						if e.IsMANE {
							vep.maneExact++
						}
					}
					for _, a := range vepAnns {
						if normalizeProteinStr(a.hgvsp) == normalizeProteinStr(expected) {
							vep.anyMatch++
							break
						}
					}
					expClass := inferConsequenceClass(expected)
					if expClass != "" {
						vep.consequTotal++
						if consequenceMatches(best.consequence, expClass) {
							vep.consequMatch++
						}
					}
				} else {
					vep.notAnnotated++
				}
			}
		}
	}

	// Write report.
	reportPath := filepath.Join(clinvarDir, "benchmark_report.md")
	if err := writeClinVarReport(reportPath, entries, vibe, snpEff, vep, hasSnpEff, hasVEP,
		vibeDur, vibeRate, cacheDur, cacheSource); err != nil {
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
	t.Logf("\nClinVar Benchmark Results (n=%d pathogenic SNVs, %d MANE Select)", len(entries), maneCount)
	t.Logf("%-10s  %s  %s  %s  %s", "Tool", "HGVSp(best)", "HGVSp(any)", "Consequence", "MANE HGVSp")
	t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s", "vibe-vep",
		pct(vibe.exactMatch, vibe.total), pct(vibe.anyMatch, vibe.total),
		pct(vibe.consequMatch, vibe.consequTotal), pct(vibe.maneExact, vibe.maneTotal))
	if hasSnpEff {
		t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s", "snpEff",
			pct(snpEff.exactMatch, snpEff.total), pct(snpEff.anyMatch, snpEff.total),
			pct(snpEff.consequMatch, snpEff.consequTotal), pct(snpEff.maneExact, snpEff.maneTotal))
	}
	if hasVEP {
		t.Logf("%-10s  %-11s  %-10s  %-11s  %-10s", "VEP",
			pct(vep.exactMatch, vep.total), pct(vep.anyMatch, vep.total),
			pct(vep.consequMatch, vep.consequTotal), pct(vep.maneExact, vep.maneTotal))
	}
}

// clinvarCounts holds benchmark counts for one annotation tool.
type clinvarCounts struct {
	total, exactMatch, anyMatch, maneTotal, maneExact int
	consequMatch, consequTotal, notAnnotated           int
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
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t)
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
// It strips the "p." prefix, normalizes "Ter" to "*", and lowercases for robustness.
// It returns the canonical form for string equality checks.
var reProteinNorm = regexp.MustCompile(`(?i)ter\b`)

func normalizeProteinStr(p string) string {
	p = strings.TrimPrefix(p, "p.")
	p = strings.TrimSpace(p)
	p = reProteinNorm.ReplaceAllString(p, "Ter")
	return p
}

// inferConsequenceClass infers the expected SO consequence class from a protein HGVS string.
func inferConsequenceClass(protein string) string {
	p := strings.TrimPrefix(protein, "p.")
	if strings.Contains(p, "fs") || strings.Contains(p, "Fs") {
		return "frameshift"
	}
	if strings.Contains(p, "del") && !strings.Contains(p, "del"+strings.ToLower(p[3:6])) {
		// p.delXxx or in-frame del — handled below
	}
	lower := strings.ToLower(p)
	switch {
	case strings.HasSuffix(lower, "ter") || strings.Contains(lower, "ter"):
		return "stop_gained"
	case strings.Contains(lower, "="):
		return "synonymous"
	case strings.Contains(lower, "del"):
		// Could be missense-adjacent or inframe
		return "" // skip ambiguous
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
	fmt.Fprintf(w, "**Ground truth**: ClinVar pathogenic/likely-pathogenic SNVs with curated protein\n")
	fmt.Fprintf(w, "HGVS annotations (from the `Name` field of `variant_summary.txt.gz`). ClinVar\n")
	fmt.Fprintf(w, "annotations are submitted by clinical labs, expert panels, and NCBI curators —\n")
	fmt.Fprintf(w, "independent of VEP, snpEff, or vibe-vep.\n\n")
	fmt.Fprintf(w, "**MANE Select**: %d / %d variants (%.0f%%) map to MANE Select transcripts.\n",
		maneCount, len(entries), 100*float64(maneCount)/float64(len(entries)))
	fmt.Fprintf(w, "For these, NM_ and ENST transcripts encode identical protein sequences,\n")
	fmt.Fprintf(w, "making cross-tool HGVSp comparison directly meaningful.\n\n")
	fmt.Fprintf(w, "**Database versions**: GENCODE v49 / Ensembl 115 (vibe-vep, VEP), snpEff GRCh38.115.\n\n")
	fmt.Fprintf(w, "## Dataset\n\n")
	fmt.Fprintf(w, "| Metric | Value |\n|--------|-------|\n")
	fmt.Fprintf(w, "| Total pathogenic SNVs | %d |\n", len(entries))
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
	fmt.Fprintf(w, "### Consequence Class Match\n\n")
	fmt.Fprintf(w, "_Inferred from protein notation: missense → `missense_variant`,_\n")
	fmt.Fprintf(w, "_Ter/\\* suffix → `stop_gained`, `fs` → `frameshift_variant`._\n\n")
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
	fmt.Fprintf(w, "| Tool | Variants | Time | Rate |\n|------|----------|------|------|\n")
	fmt.Fprintf(w, "| vibe-vep | %d | %.1fs | %.0f v/s |\n", len(entries), vibeDur.Seconds(), vibeRate)
	fmt.Fprintf(w, "_Cache load: %.1fs from %s_\n\n", cacheDur.Seconds(), cacheSource)
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
	fmt.Fprintf(w, "providing the most rigorous comparison of actual prediction accuracy.\n")

	return w.Flush()
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
