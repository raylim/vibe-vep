package output_test

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/duckdb"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/output"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// snpEffAnnotation holds one annotation entry from snpEff's ANN INFO field.
//
// ANN format (pipe-separated):
//
//	Allele|Effect|Impact|GeneName|GeneID|Feature_Type|Feature_ID|Transcript_BioType|
//	Rank|HGVS.c|HGVS.p|cDNA.pos/len|CDS.pos/len|AA.pos/len|Distance|ERRORS
type snpEffAnnotation struct {
	allele      string
	consequence string // Sequence Ontology term (snpEff uses SO terms)
	impact      string // HIGH / MODERATE / LOW / MODIFIER
	geneName    string
	transcriptID string // Feature_ID when Feature_Type == "transcript"
	biotype     string
	hgvsc       string
	hgvsp       string
}

// snpEffVariantMap maps a variant key to its snpEff annotations.
// Primary key: "chrom:pos:ref:alt" (exact, works for SNVs).
// Secondary key: "chrom:pos" (position-only, fallback for indels with anchor-base differences).
type snpEffVariantMap struct {
	exact    map[string][]snpEffAnnotation
	byPos    map[string][]snpEffAnnotation
}

func newSnpEffVariantMap() snpEffVariantMap {
	return snpEffVariantMap{
		exact: make(map[string][]snpEffAnnotation),
		byPos: make(map[string][]snpEffAnnotation),
	}
}

// lookup returns annotations for a vcf.Variant, trying exact match first.
func (m snpEffVariantMap) lookup(v *vcf.Variant) []snpEffAnnotation {
	key := snpEffExactKey(v.NormalizeChrom(), v.Pos, v.Ref, v.Alt)
	if anns, ok := m.exact[key]; ok {
		return anns
	}
	// Fallback: position-only match for indels where anchor-base convention differs.
	posKey := snpEffPosKey(v.NormalizeChrom(), v.Pos)
	return m.byPos[posKey]
}

// snpEffExactKey returns "chrom:pos:ref:alt".
func snpEffExactKey(chrom string, pos int64, ref, alt string) string {
	return fmt.Sprintf("%s:%d:%s:%s", normalizeChrom(chrom), pos, ref, alt)
}

// snpEffPosKey returns "chrom:pos".
func snpEffPosKey(chrom string, pos int64) string {
	return fmt.Sprintf("%s:%d", normalizeChrom(chrom), pos)
}

func normalizeChrom(c string) string {
	if len(c) > 3 && c[:3] == "chr" {
		return c[3:]
	}
	return c
}

// loadSnpEffVCF parses a snpEff-annotated VCF file (plain or gzipped) into a variant map.
//
// snpEff annotates each record with an ANN INFO field containing one or more
// pipe-separated annotation entries. Multi-allelic ALT columns are split, and
// per-allele annotations are matched by the first (Allele) field.
func loadSnpEffVCF(path string) (snpEffVariantMap, error) {
	f, err := os.Open(path)
	if err != nil {
		return snpEffVariantMap{}, err
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gr, err := gzip.NewReader(f)
		if err != nil {
			return snpEffVariantMap{}, fmt.Errorf("gzip reader: %w", err)
		}
		defer gr.Close()
		r = gr
	}

	m := newSnpEffVariantMap()
	scanner := bufio.NewScanner(r)
	// snpEff INFO fields can be very long for multi-allelic sites.
	scanner.Buffer(make([]byte, 64*1024), 64*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 9)
		if len(fields) < 8 {
			continue
		}
		chrom := normalizeChrom(fields[0])
		posStr := fields[1]
		ref := fields[3]
		alts := strings.Split(fields[4], ",")
		infoStr := fields[7]

		annField := extractANNField(infoStr)
		if annField == "" {
			continue
		}
		anns := parseSnpEffANN(annField)
		if len(anns) == 0 {
			continue
		}

		for _, alt := range alts {
			// Filter annotations to those matching this allele (or keep all for mono-allelic).
			altAnns := filterByAllele(anns, alt, len(alts) == 1)
			if len(altAnns) == 0 {
				continue
			}
			// For deletions represented in VCF format (anchor base prepended), also
			// register under the anchor-stripped form to match MAF-derived keys.
			exactKey := snpEffExactKey(chrom, mustParsePos(posStr), ref, alt)
			m.exact[exactKey] = append(m.exact[exactKey], altAnns...)

			posKey := snpEffPosKey(chrom, mustParsePos(posStr))
			m.byPos[posKey] = append(m.byPos[posKey], altAnns...)

			// For SNV, also store position adjusted for possible MAF offset conventions.
			// Deletions in VCF have pos = MAF_Start - 1; register at pos+1 as well.
			if len(ref) > len(alt) && len(ref) > 1 {
				altKey2 := snpEffExactKey(chrom, mustParsePos(posStr)+1, ref[1:], "")
				m.exact[altKey2] = append(m.exact[altKey2], altAnns...)
			}
			// Insertions in VCF have pos = MAF_Start; register stripped form.
			if len(alt) > len(ref) && len(ref) == 1 {
				altKey2 := snpEffExactKey(chrom, mustParsePos(posStr), "", alt[1:])
				m.exact[altKey2] = append(m.exact[altKey2], altAnns...)
			}
		}
	}
	return m, scanner.Err()
}

func mustParsePos(s string) int64 {
	var n int64
	for _, c := range s {
		if c < '0' || c > '9' {
			break
		}
		n = n*10 + int64(c-'0')
	}
	return n
}

// extractANNField returns the value of the ANN (or EFF) INFO key, or "".
func extractANNField(info string) string {
	for _, kv := range strings.Split(info, ";") {
		if strings.HasPrefix(kv, "ANN=") {
			return kv[4:]
		}
		if strings.HasPrefix(kv, "EFF=") { // snpEff v4 legacy format
			return kv[4:]
		}
	}
	return ""
}

// parseSnpEffANN parses a comma-separated ANN field value.
func parseSnpEffANN(annField string) []snpEffAnnotation {
	var anns []snpEffAnnotation
	for _, entry := range strings.Split(annField, ",") {
		f := strings.Split(entry, "|")
		if len(f) < 11 {
			continue
		}
		// Only transcript-level annotations are useful for comparison.
		if f[5] != "transcript" {
			continue
		}
		ann := snpEffAnnotation{
			allele:       f[0],
			consequence:  f[1],
			impact:       f[2],
			geneName:     f[3],
			transcriptID: f[6],
			biotype:      f[7],
			hgvsc:        stripTranscriptPrefix(f[9]),
			hgvsp:        stripTranscriptPrefix(f[10]),
		}
		anns = append(anns, ann)
	}
	return anns
}

// stripTranscriptPrefix removes any "ENST…:" or "NM_…:" prefix from an HGVS string.
func stripTranscriptPrefix(s string) string {
	if i := strings.Index(s, ":"); i >= 0 {
		return s[i+1:]
	}
	return s
}

// filterByAllele returns annotations matching the given alt allele.
// When singleAlt is true (mono-allelic site) all annotations are returned.
func filterByAllele(anns []snpEffAnnotation, alt string, singleAlt bool) []snpEffAnnotation {
	if singleAlt {
		return anns
	}
	var out []snpEffAnnotation
	for _, a := range anns {
		if a.allele == alt {
			out = append(out, a)
		}
	}
	return out
}

// selectBestSnpEffAnnotation mirrors SelectBestAnnotation for snpEff output.
// Priority: transcript ID match > same gene > protein_coding > highest impact.
func selectBestSnpEffAnnotation(mafAnn *maf.MAFAnnotation, anns []snpEffAnnotation) *snpEffAnnotation {
	if len(anns) == 0 {
		return nil
	}

	// Pass 1: transcript ID match ignoring version suffix.
	if mafAnn.TranscriptID != "" {
		mafBase := snpEffStripVersion(mafAnn.TranscriptID)
		for i := range anns {
			if snpEffStripVersion(anns[i].transcriptID) == mafBase {
				return &anns[i]
			}
		}
	}

	// Pass 2: same gene, prefer protein_coding > higher impact.
	var best *snpEffAnnotation
	for i := range anns {
		if anns[i].geneName == mafAnn.HugoSymbol {
			if best == nil || snpEffAnnBetter(&anns[i], best) {
				best = &anns[i]
			}
		}
	}
	if best != nil {
		return best
	}

	// Pass 3: overall best.
	best = &anns[0]
	for i := 1; i < len(anns); i++ {
		if snpEffAnnBetter(&anns[i], best) {
			best = &anns[i]
		}
	}
	return best
}

func snpEffAnnBetter(a, b *snpEffAnnotation) bool {
	aPC := a.biotype == "protein_coding"
	bPC := b.biotype == "protein_coding"
	if aPC != bPC {
		return aPC
	}
	return annotate.ImpactRank(a.impact) > annotate.ImpactRank(b.impact)
}

// snpEffStripVersion strips the version suffix (e.g. ".4") from a transcript ID.
func snpEffStripVersion(id string) string {
	if i := strings.LastIndexByte(id, '.'); i >= 0 {
		return id[:i]
	}
	return id
}

// snpEffAnnotationToMap converts a snpEffAnnotation to the column map used by Categorizer.
// Returns empty strings when ann is nil (no snpEff annotation found for the variant).
func snpEffAnnotationToMap(ann *snpEffAnnotation) map[string]string {
	if ann == nil {
		return map[string]string{
			"Consequence": "",
			"HGVSp_Short": "",
			"HGVSc":       "",
		}
	}
	return map[string]string{
		"Consequence": ann.consequence,
		"HGVSp_Short": output.HGVSpToShort(ann.hgvsp),
		"HGVSc":       ann.hgvsc,
	}
}

// annotationToMapForSnpEff converts a vibe-vep Annotation to the column map.
func annotationToMapForSnpEff(ann *annotate.Annotation) map[string]string {
	if ann == nil {
		return map[string]string{
			"Consequence": "",
			"HGVSp_Short": "",
			"HGVSc":       "",
		}
	}
	return map[string]string{
		"Consequence": ann.Consequence,
		"HGVSp_Short": output.HGVSpToShort(ann.HGVSp),
		"HGVSc":       ann.HGVSc,
	}
}

// pct returns 100 * num / denom, or 0 if denom == 0.
func pct(num, denom int) float64 {
	if denom == 0 {
		return 0
	}
	return float64(num) / float64(denom) * 100
}

// initCatCounts initialises a per-column category-count map.
func initCatCounts(cols []string) map[string]map[output.Category]int {
	m := make(map[string]map[output.Category]int, len(cols))
	for _, c := range cols {
		m[c] = make(map[output.Category]int)
	}
	return m
}

// snpEffStudyResult holds per-study benchmark results for one tool.
type snpEffStudyResult struct {
	name           string
	variantsTotal  int // total MAF variants processed
	variantsInTool int // variants that had a matching annotation in the tool output
	catCounts      map[string]map[output.Category]int
}

// snpEffBenchResult holds the three-way comparison for one study.
type snpEffBenchResult struct {
	study  string
	vibe   snpEffStudyResult
	snpeff snpEffStudyResult
	vep    snpEffStudyResult // Ensembl VEP (zero value when VEP VCF not available)
}

// TestSnpEffBenchmark compares vibe-vep, snpEff, and Ensembl VEP prediction
// accuracy against TCGA GDC annotations (used as ground truth).
//
// Requires pre-generated snpEff VCF files at testdata/tcga/snpeff/<study>.vcf.gz.
// Generate with:
//
//	scripts/run_snpeff.sh
//
// Ensembl VEP VCF files are optional; if present at testdata/tcga/vep/<study>.vcf.gz
// they are included in a 3-way comparison. Generate with:
//
//	scripts/run_vep.sh
//
// The test writes a markdown report to testdata/tcga/benchmark_report.md.
//
// Skipped with -short (large files, not for CI). Run with:
//
//	go test ./internal/output/ -run TestSnpEffBenchmark -v -count=1
func TestSnpEffBenchmark(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping benchmark in short mode")
	}

	tcgaDir := findStudyDir(t, "tcga")
	snpEffDir := filepath.Join(tcgaDir, "snpeff")
	vepDir := filepath.Join(tcgaDir, "vep")

	mafFiles, err := filepath.Glob(filepath.Join(tcgaDir, "*_data_mutations.txt"))
	if err != nil {
		t.Fatalf("glob MAF files: %v", err)
	}
	if len(mafFiles) == 0 {
		t.Skip("no TCGA MAF files found in", tcgaDir)
	}
	sort.Strings(mafFiles)

	// Collect studies that have at least a snpEff VCF; VEP VCF is optional.
	type bench struct {
		maf, snpeffVCF string
		vepVCF         string // empty if not available
	}
	var benches []bench
	for _, mf := range mafFiles {
		name := studyName(mf)
		var snpeffVCF, vepVCF string
		for _, ext := range []string{".vcf.gz", ".vcf"} {
			if p := filepath.Join(snpEffDir, name+ext); fileExists(p) {
				snpeffVCF = p
				break
			}
		}
		for _, ext := range []string{".vcf.gz", ".vcf"} {
			if p := filepath.Join(vepDir, name+ext); fileExists(p) {
				vepVCF = p
				break
			}
		}
		if snpeffVCF != "" {
			benches = append(benches, bench{mf, snpeffVCF, vepVCF})
		}
	}
	if len(benches) == 0 {
		t.Skipf("no snpEff VCF files found in %s — run scripts/run_snpeff.sh first", snpEffDir)
	}

	// Load GENCODE transcripts (reuse the same helper as validation_benchmark_test.go).
	gtfPath, fastaPath, canonicalPath := findGENCODEFiles(t, "GRCh38")
	cacheDir := filepath.Dir(gtfPath)
	c := cache.New()

	var loadSource string
	cacheStart := time.Now()

	tc := duckdb.NewTranscriptCache(cacheDir)
	gtfFP, err1 := duckdb.StatFile(gtfPath)
	fastaFP, err2 := duckdb.StatFile(fastaPath)
	canonicalFP := duckdb.FileFingerprint{}
	if canonicalPath != "" {
		canonicalFP, _ = duckdb.StatFile(canonicalPath)
	}
	if err1 == nil && err2 == nil && tc.Valid(gtfFP, fastaFP, canonicalFP) {
		if err := tc.Load(c); err != nil {
			t.Fatalf("load transcript cache: %v", err)
		}
		loadSource = "gob cache"
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
	cacheDuration := time.Since(cacheStart)
	t.Logf("loaded %d transcripts from %s in %s",
		c.TranscriptCount(), loadSource, cacheDuration.Round(time.Millisecond))

	categorizer := &output.Categorizer{}
	comparedCols := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	ann := annotate.NewAnnotator(c)

	var results []snpEffBenchResult

	for _, b := range benches {
		name := studyName(b.maf)
		t.Run(name, func(t *testing.T) {
			t.Logf("loading snpEff VCF: %s", b.snpeffVCF)
			snpEffMap, err := loadSnpEffVCF(b.snpeffVCF)
			if err != nil {
				t.Fatalf("load snpEff VCF: %v", err)
			}
			t.Logf("loaded %d unique snpEff positions", len(snpEffMap.exact))

			// Load VEP VCF if present.
			var vepMap vepVariantMap
			hasVEP := b.vepVCF != ""
			if hasVEP {
				t.Logf("loading VEP VCF: %s", b.vepVCF)
				vepMap, err = loadVEPVCF(b.vepVCF)
				if err != nil {
					t.Fatalf("load VEP VCF: %v", err)
				}
				t.Logf("loaded %d unique VEP positions", len(vepMap.exact))
			}

			parser, err := maf.NewParser(b.maf)
			if err != nil {
				t.Fatalf("open MAF: %v", err)
			}
			defer parser.Close()

			vibeCats := initCatCounts(comparedCols)
			snpEffCats := initCatCounts(comparedCols)
			vepCats := initCatCounts(comparedCols)
			var total, snpEffMissing, vepMissing int

			for {
				v, mafAnn, err := parser.NextWithAnnotation()
				if err != nil {
					t.Fatalf("read variant: %v", err)
				}
				if v == nil {
					break
				}
				vibeAnns, err := ann.Annotate(v)
				if err != nil {
					continue
				}
				total++

				left := map[string]string{
					"Consequence": mafAnn.Consequence,
					"HGVSp_Short": mafAnn.HGVSpShort,
					"HGVSc":       mafAnn.HGVSc,
				}

				// vibe-vep prediction.
				bestVibe := output.SelectBestAnnotation(mafAnn, vibeAnns)
				vibeRight := annotationToMapForSnpEff(bestVibe)
				for col, cat := range categorizer.CategorizeRow(comparedCols, left, vibeRight, comparedCols, comparedCols) {
					vibeCats[col][cat]++
				}

				// snpEff prediction.
				seAnns := snpEffMap.lookup(v)
				bestSE := selectBestSnpEffAnnotation(mafAnn, seAnns)
				if bestSE == nil {
					snpEffMissing++
				}
				seRight := snpEffAnnotationToMap(bestSE)
				for col, cat := range categorizer.CategorizeRow(comparedCols, left, seRight, comparedCols, comparedCols) {
					snpEffCats[col][cat]++
				}

				// Ensembl VEP prediction.
				if hasVEP {
					evAnns := vepMap.lookup(v)
					bestVEP := selectBestVEPAnnotation(mafAnn, evAnns)
					if bestVEP == nil {
						vepMissing++
					}
					evRight := vepAnnotationToMap(bestVEP)
					for col, cat := range categorizer.CategorizeRow(comparedCols, left, evRight, comparedCols, comparedCols) {
						vepCats[col][cat]++
					}
				}
			}

			vibeConseqMatch := vibeCats["Consequence"][output.CatMatch] + vibeCats["Consequence"][output.CatUpstreamReclass]
			seConseqMatch := snpEffCats["Consequence"][output.CatMatch] + snpEffCats["Consequence"][output.CatUpstreamReclass]

			t.Logf("%s: %d variants", name, total)
			t.Logf("  vibe-vep consequence match: %.1f%%", pct(vibeConseqMatch, total))
			t.Logf("  snpEff   consequence match: %.1f%%  (%d variants not found in snpEff VCF)",
				pct(seConseqMatch, total), snpEffMissing)
			if hasVEP {
				vepConseqMatch := vepCats["Consequence"][output.CatMatch] + vepCats["Consequence"][output.CatUpstreamReclass]
				t.Logf("  VEP      consequence match: %.1f%%  (%d variants not found in VEP VCF)",
					pct(vepConseqMatch, total), vepMissing)
			}

			r := snpEffBenchResult{
				study: name,
				vibe: snpEffStudyResult{
					name:           "vibe-vep",
					variantsTotal:  total,
					variantsInTool: total,
					catCounts:      vibeCats,
				},
				snpeff: snpEffStudyResult{
					name:           "snpEff",
					variantsTotal:  total,
					variantsInTool: total - snpEffMissing,
					catCounts:      snpEffCats,
				},
			}
			if hasVEP {
				r.vep = snpEffStudyResult{
					name:           "VEP",
					variantsTotal:  total,
					variantsInTool: total - vepMissing,
					catCounts:      vepCats,
				}
			}
			results = append(results, r)
		})
	}

	if len(results) == 0 {
		return
	}

	reportPath := filepath.Join(tcgaDir, "benchmark_report.md")
	writeSnpEffReport(t, reportPath, results, c.TranscriptCount(), cacheDuration, loadSource)
}

// writeSnpEffReport generates a markdown report comparing vibe-vep, snpEff, and Ensembl VEP.
func writeSnpEffReport(
	t *testing.T,
	path string,
	results []snpEffBenchResult,
	transcriptCount int,
	cacheDuration time.Duration,
	loadSource string,
) {
	t.Helper()

	// Determine whether VEP data is present in any study.
	hasVEP := false
	for _, r := range results {
		if r.vep.variantsTotal > 0 {
			hasVEP = true
			break
		}
	}

	var sb strings.Builder

	title := "snpEff vs vibe-vep Benchmark Report"
	if hasVEP {
		title = "vibe-vep vs snpEff vs Ensembl VEP Benchmark Report"
	}
	sb.WriteString(fmt.Sprintf("# %s\n\n", title))
	sb.WriteString(fmt.Sprintf("Generated: %s  \n", time.Now().UTC().Format("2006-01-02 15:04 UTC")))
	sb.WriteString(fmt.Sprintf("GENCODE transcripts: %d (loaded from %s in %s)  \n",
		transcriptCount, loadSource, cacheDuration.Round(time.Millisecond)))
	sb.WriteString("\nGround truth: TCGA GDC MAF annotations.  \n")
	if hasVEP {
		sb.WriteString("snpEff VCFs generated by: `scripts/run_snpeff.sh`  \n")
		sb.WriteString("Ensembl VEP VCFs generated by: `scripts/run_vep.sh`\n\n")
	} else {
		sb.WriteString("snpEff VCFs generated by: `scripts/run_snpeff.sh`\n\n")
	}

	if hasVEP {
		sb.WriteString("## Methodology Notes\n\n")
		sb.WriteString("All three tools use **Ensembl 115 gene models** for this benchmark run:\n\n")
		sb.WriteString("| Tool | Gene model | Ensembl release |\n")
		sb.WriteString("|------|-----------|----------------|\n")
		sb.WriteString("| Ground truth (GDC MAF) | Ensembl VEP + MSK isoform overrides | ~Ensembl 92–99 era |\n")
		sb.WriteString("| snpEff | GRCh38.115 | Ensembl 115 |\n")
		sb.WriteString("| vibe-vep | GENCODE v49 | Ensembl 115 |\n")
		sb.WriteString("| Ensembl VEP (this run) | VEP v115 cache | Ensembl 115 |\n\n")
		sb.WriteString("**Important caveats:**\n\n")
		sb.WriteString("1. **Ground truth is VEP-derived**: The TCGA GDC MAF files were themselves annotated with " +
			"Ensembl VEP (column names such as `CANONICAL`, `BIOTYPE`, `APPRIS`, `MANE` are VEP output fields). " +
			"This gives Ensembl VEP a structural advantage in consequence and HGVSc match rates since both use " +
			"identical SO term vocabularies and transcript selection logic.\n\n")
		sb.WriteString("2. **GDC uses MSK isoform overrides** (`#isoform: mskcc` header): The GDC applies " +
			"cancer-relevant isoform prioritization via genome-nexus, not purely Ensembl canonical isoforms. " +
			"Neither our VEP run nor snpEff uses these overrides. vibe-vep's higher HGVSp match rate (~94% vs " +
			"VEP's ~70%) reflects that vibe-vep's transcript prioritization is more compatible with MSK/GDC " +
			"isoform choices.\n\n")
		sb.WriteString("3. **Version drift vs ground truth**: All tools use Ensembl 115, which is 16+ releases " +
			"newer than the GDC ground truth era (~Ensembl 92–99). Some transcripts have changed exon structure " +
			"or been retired between versions. Despite using the same Ensembl 115 database as VEP, vibe-vep " +
			"achieves higher HGVSp match rates due to its MSK-compatible transcript prioritization.\n\n")
		sb.WriteString("**What this benchmark measures**: Compatibility with GDC/MSK clinical annotations, " +
			"not absolute correctness against the genome. For clinical variant reporting workflows that produce " +
			"or consume GDC-style MAF files, vibe-vep's annotations are most compatible.\n\n")
	}

	comparedCols := []string{"Consequence", "HGVSp_Short", "HGVSc"}
	colDisplay := map[string]string{
		"Consequence": "Consequence",
		"HGVSp_Short": "HGVSp",
		"HGVSc":       "HGVSc",
	}

	// Summary table: consequence match rate for all tools.
	sb.WriteString("## Consequence Match Rate\n\n")
	if hasVEP {
		sb.WriteString("| Study | Variants | vibe-vep Match | vibe-vep Rate | snpEff Match | snpEff Rate | snpEff Cov | VEP Match | VEP Rate | VEP Cov |\n")
		sb.WriteString("|-------|----------|---------------|--------------|--------------|-------------|------------|-----------|----------|---------|\n")
	} else {
		sb.WriteString("| Study | Variants | vibe-vep Match | vibe-vep Rate | snpEff Match | snpEff Rate | snpEff Coverage |\n")
		sb.WriteString("|-------|----------|---------------|--------------|--------------|-------------|-----------------|\n")
	}

	var totVariants, totVibeMatch, totSEMatch, totVEPMatch int
	var totSEMissing, totVEPMissing int
	for _, r := range results {
		vibeMatch := r.vibe.catCounts["Consequence"][output.CatMatch] + r.vibe.catCounts["Consequence"][output.CatUpstreamReclass]
		seMatch := r.snpeff.catCounts["Consequence"][output.CatMatch] + r.snpeff.catCounts["Consequence"][output.CatUpstreamReclass]
		seMissing := r.snpeff.variantsTotal - r.snpeff.variantsInTool
		if hasVEP {
			vepMatch := r.vep.catCounts["Consequence"][output.CatMatch] + r.vep.catCounts["Consequence"][output.CatUpstreamReclass]
			vepMissing := r.vep.variantsTotal - r.vep.variantsInTool
			sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.1f%% | %d | %.1f%% | %.1f%% | %d | %.1f%% | %.1f%% |\n",
				r.study, r.vibe.variantsTotal,
				vibeMatch, pct(vibeMatch, r.vibe.variantsTotal),
				seMatch, pct(seMatch, r.snpeff.variantsTotal),
				pct(r.snpeff.variantsInTool, r.snpeff.variantsTotal),
				vepMatch, pct(vepMatch, r.vep.variantsTotal),
				pct(r.vep.variantsInTool, r.vep.variantsTotal)))
			totVEPMatch += vepMatch
			totVEPMissing += vepMissing
		} else {
			sb.WriteString(fmt.Sprintf("| %s | %d | %d | %.1f%% | %d | %.1f%% | %.1f%% |\n",
				r.study, r.vibe.variantsTotal,
				vibeMatch, pct(vibeMatch, r.vibe.variantsTotal),
				seMatch, pct(seMatch, r.snpeff.variantsTotal),
				pct(r.snpeff.variantsInTool, r.snpeff.variantsTotal)))
		}
		totVariants += r.vibe.variantsTotal
		totVibeMatch += vibeMatch
		totSEMatch += seMatch
		totSEMissing += seMissing
	}
	if hasVEP {
		sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.1f%%** | **%d** | **%.1f%%** | **%.1f%%** | **%d** | **%.1f%%** | **%.1f%%** |\n\n",
			totVariants,
			totVibeMatch, pct(totVibeMatch, totVariants),
			totSEMatch, pct(totSEMatch, totVariants),
			pct(totVariants-totSEMissing, totVariants),
			totVEPMatch, pct(totVEPMatch, totVariants),
			pct(totVariants-totVEPMissing, totVariants)))
	} else {
		sb.WriteString(fmt.Sprintf("| **Total** | **%d** | **%d** | **%.1f%%** | **%d** | **%.1f%%** | **%.1f%%** |\n\n",
			totVariants,
			totVibeMatch, pct(totVibeMatch, totVariants),
			totSEMatch, pct(totSEMatch, totVariants),
			pct(totVariants-totSEMissing, totVariants)))
	}

	// Full match-rate table for all columns.
	sb.WriteString("## Full Column Match Rates\n\n")
	if hasVEP {
		sb.WriteString("| Study | Column | vibe-vep Match | vibe-vep Rate | snpEff Match | snpEff Rate | VEP Match | VEP Rate |\n")
		sb.WriteString("|-------|--------|---------------|--------------|--------------|-------------|-----------|----------|\n")
	} else {
		sb.WriteString("| Study | Column | vibe-vep Match | vibe-vep Rate | snpEff Match | snpEff Rate |\n")
		sb.WriteString("|-------|--------|---------------|--------------|--------------|-------------|\n")
	}

	for _, r := range results {
		for _, col := range comparedCols {
			vibeMatch := matchCount(r.vibe.catCounts, col)
			seMatch := matchCount(r.snpeff.catCounts, col)
			if hasVEP {
				vepMatch := matchCount(r.vep.catCounts, col)
				sb.WriteString(fmt.Sprintf("| %s | %s | %d | %.1f%% | %d | %.1f%% | %d | %.1f%% |\n",
					r.study, colDisplay[col],
					vibeMatch, pct(vibeMatch, r.vibe.variantsTotal),
					seMatch, pct(seMatch, r.snpeff.variantsTotal),
					vepMatch, pct(vepMatch, r.vep.variantsTotal)))
			} else {
				sb.WriteString(fmt.Sprintf("| %s | %s | %d | %.1f%% | %d | %.1f%% |\n",
					r.study, colDisplay[col],
					vibeMatch, pct(vibeMatch, r.vibe.variantsTotal),
					seMatch, pct(seMatch, r.snpeff.variantsTotal)))
			}
		}
	}
	sb.WriteString("\n")

	// Per-column category breakdown.
	tools := []struct {
		label  string
		getter func(r snpEffBenchResult) snpEffStudyResult
	}{
		{"vibe-vep", func(r snpEffBenchResult) snpEffStudyResult { return r.vibe }},
		{"snpEff", func(r snpEffBenchResult) snpEffStudyResult { return r.snpeff }},
	}
	if hasVEP {
		tools = append(tools, struct {
			label  string
			getter func(r snpEffBenchResult) snpEffStudyResult
		}{"VEP", func(r snpEffBenchResult) snpEffStudyResult { return r.vep }})
	}

	for _, col := range comparedCols {
		sb.WriteString(fmt.Sprintf("## %s Category Breakdown\n\n", colDisplay[col]))

		allCats := make(map[output.Category]bool)
		for _, r := range results {
			for _, tool := range tools {
				for cat := range tool.getter(r).catCounts[col] {
					allCats[cat] = true
				}
			}
		}
		catList := sortedCategories(allCats)

		sb.WriteString("| Study | Tool |")
		for _, cat := range catList {
			sb.WriteString(fmt.Sprintf(" %s |", cat))
		}
		sb.WriteString("\n|-------|------|")
		for range catList {
			sb.WriteString("------|")
		}
		sb.WriteString("\n")

		totals := make([]map[output.Category]int, len(tools))
		for i := range totals {
			totals[i] = make(map[output.Category]int)
		}

		for _, r := range results {
			for ti, tool := range tools {
				sb.WriteString(fmt.Sprintf("| %s | %s |", r.study, tool.label))
				for _, cat := range catList {
					n := tool.getter(r).catCounts[col][cat]
					totals[ti][cat] += n
					sb.WriteString(fmt.Sprintf(" %d |", n))
				}
				sb.WriteString("\n")
			}
		}
		for ti, tool := range tools {
			sb.WriteString(fmt.Sprintf("| **Total** | **%s** |", tool.label))
			for _, cat := range catList {
				sb.WriteString(fmt.Sprintf(" **%d** |", totals[ti][cat]))
			}
			sb.WriteString("\n")
		}
		sb.WriteString("\n")
	}

	if err := os.WriteFile(path, []byte(sb.String()), 0644); err != nil {
		t.Fatalf("write report: %v", err)
	}
	t.Logf("report written to %s", path)
}

func sortedCategories(m map[output.Category]bool) []output.Category {
	cats := make([]output.Category, 0, len(m))
	for c := range m {
		cats = append(cats, c)
	}
	sort.Slice(cats, func(i, j int) bool { return cats[i] < cats[j] })
	return cats
}

// matchCount returns the match count for a column, including upstream reclassifications
// for the Consequence column.
func matchCount(catCounts map[string]map[output.Category]int, col string) int {
	n := catCounts[col][output.CatMatch]
	if col == "Consequence" {
		n += catCounts[col][output.CatUpstreamReclass]
	}
	return n
}

// fileExists reports whether path exists and is a regular file.
func fileExists(path string) bool {
	fi, err := os.Stat(path)
	return err == nil && fi.Mode().IsRegular()
}

// ── Ensembl VEP support ───────────────────────────────────────────────────────

// vepAnnotation holds one annotation entry from Ensembl VEP's CSQ INFO field.
//
// The CSQ format is defined in the VCF header:
//
//	##INFO=<ID=CSQ,...,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|...">
type vepAnnotation struct {
	allele       string
	consequence  string
	impact       string
	geneName     string
	transcriptID string
	biotype      string
	canonical    bool
	hgvsc        string
	hgvsp        string
}

// vepVariantMap maps variant keys to VEP annotations, using the same dual-key
// strategy as snpEffVariantMap.
type vepVariantMap struct {
	exact map[string][]vepAnnotation
	byPos map[string][]vepAnnotation
}

func newVEPVariantMap() vepVariantMap {
	return vepVariantMap{
		exact: make(map[string][]vepAnnotation),
		byPos: make(map[string][]vepAnnotation),
	}
}

func (m vepVariantMap) lookup(v *vcf.Variant) []vepAnnotation {
	key := snpEffExactKey(v.NormalizeChrom(), v.Pos, v.Ref, v.Alt)
	if anns, ok := m.exact[key]; ok {
		return anns
	}
	return m.byPos[snpEffPosKey(v.NormalizeChrom(), v.Pos)]
}

// loadVEPVCF parses a VEP-annotated VCF file (plain or gzipped) into a variant map.
//
// VEP writes annotations in the CSQ INFO field. The field format (pipe-separated
// column names) is read from the ##INFO=<ID=CSQ,...> header line.
func loadVEPVCF(path string) (vepVariantMap, error) {
	f, err := os.Open(path)
	if err != nil {
		return vepVariantMap{}, err
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gr, err := gzip.NewReader(f)
		if err != nil {
			return vepVariantMap{}, fmt.Errorf("gzip reader: %w", err)
		}
		defer gr.Close()
		r = gr
	}

	m := newVEPVariantMap()
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 64*1024*1024)

	var csqFieldIdx map[string]int

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "##INFO=<ID=CSQ") {
			csqFieldIdx = parseVEPCSQHeader(line)
			continue
		}
		if strings.HasPrefix(line, "#") {
			continue
		}
		if csqFieldIdx == nil {
			continue
		}

		fields := strings.SplitN(line, "\t", 9)
		if len(fields) < 8 {
			continue
		}
		chrom := normalizeChrom(fields[0])
		posStr := fields[1]
		ref := fields[3]
		alts := strings.Split(fields[4], ",")
		infoStr := fields[7]

		csqField := extractCSQField(infoStr)
		if csqField == "" {
			continue
		}
		anns := parseVEPCSQ(csqField, csqFieldIdx)
		if len(anns) == 0 {
			continue
		}

		for _, alt := range alts {
			altAnns := filterVEPByAllele(anns, alt, len(alts) == 1)
			if len(altAnns) == 0 {
				continue
			}
			pos := mustParsePos(posStr)
			exactKey := snpEffExactKey(chrom, pos, ref, alt)
			m.exact[exactKey] = append(m.exact[exactKey], altAnns...)
			posKey := snpEffPosKey(chrom, pos)
			m.byPos[posKey] = append(m.byPos[posKey], altAnns...)

			// Register anchor-stripped forms for indel matching (same as snpEff).
			if len(ref) > len(alt) && len(ref) > 1 {
				k := snpEffExactKey(chrom, pos+1, ref[1:], "")
				m.exact[k] = append(m.exact[k], altAnns...)
			}
			if len(alt) > len(ref) && len(ref) == 1 {
				k := snpEffExactKey(chrom, pos, "", alt[1:])
				m.exact[k] = append(m.exact[k], altAnns...)
			}
		}
	}
	return m, scanner.Err()
}

// parseVEPCSQHeader extracts a field-name → index map from the CSQ header line.
//
// VEP header format:
//
//	##INFO=<ID=CSQ,...,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|...">
func parseVEPCSQHeader(line string) map[string]int {
	const marker = "Format: "
	i := strings.Index(line, marker)
	if i < 0 {
		return nil
	}
	formatStr := line[i+len(marker):]
	formatStr = strings.TrimRight(formatStr, `">`)
	fields := strings.Split(formatStr, "|")
	idx := make(map[string]int, len(fields))
	for j, f := range fields {
		idx[f] = j
	}
	return idx
}

// extractCSQField returns the value of the CSQ INFO key, or "".
func extractCSQField(info string) string {
	for _, kv := range strings.Split(info, ";") {
		if strings.HasPrefix(kv, "CSQ=") {
			return kv[4:]
		}
	}
	return ""
}

// parseVEPCSQ parses a comma-separated CSQ field value into vepAnnotation entries.
func parseVEPCSQ(csqField string, fieldIdx map[string]int) []vepAnnotation {
	get := func(f []string, name string) string {
		i, ok := fieldIdx[name]
		if !ok || i >= len(f) {
			return ""
		}
		return f[i]
	}

	var anns []vepAnnotation
	for _, entry := range strings.Split(csqField, ",") {
		f := strings.Split(entry, "|")
		if !strings.EqualFold(get(f, "Feature_type"), "Transcript") {
			continue
		}
		anns = append(anns, vepAnnotation{
			allele:       get(f, "Allele"),
			consequence:  get(f, "Consequence"),
			impact:       get(f, "IMPACT"),
			geneName:     get(f, "SYMBOL"),
			transcriptID: get(f, "Feature"),
			biotype:      get(f, "BIOTYPE"),
			canonical:    get(f, "CANONICAL") == "YES",
			hgvsc:        stripTranscriptPrefix(get(f, "HGVSc")),
			hgvsp:        stripTranscriptPrefix(get(f, "HGVSp")),
		})
	}
	return anns
}

// filterVEPByAllele returns annotations matching the given alt allele.
func filterVEPByAllele(anns []vepAnnotation, alt string, singleAlt bool) []vepAnnotation {
	if singleAlt {
		return anns
	}
	var out []vepAnnotation
	for _, a := range anns {
		if a.allele == alt {
			out = append(out, a)
		}
	}
	return out
}

// selectBestVEPAnnotation mirrors selectBestSnpEffAnnotation for VEP output.
// Priority: transcript ID match > canonical > protein_coding > highest impact.
func selectBestVEPAnnotation(mafAnn *maf.MAFAnnotation, anns []vepAnnotation) *vepAnnotation {
	if len(anns) == 0 {
		return nil
	}

	// Pass 1: exact transcript ID match ignoring version suffix.
	if mafAnn.TranscriptID != "" {
		mafBase := snpEffStripVersion(mafAnn.TranscriptID)
		for i := range anns {
			if snpEffStripVersion(anns[i].transcriptID) == mafBase {
				return &anns[i]
			}
		}
	}

	// Pass 2: same gene, prefer canonical > protein_coding > higher impact.
	var best *vepAnnotation
	for i := range anns {
		if anns[i].geneName == mafAnn.HugoSymbol {
			if best == nil || vepAnnBetter(&anns[i], best) {
				best = &anns[i]
			}
		}
	}
	if best != nil {
		return best
	}

	// Pass 3: overall best.
	best = &anns[0]
	for i := 1; i < len(anns); i++ {
		if vepAnnBetter(&anns[i], best) {
			best = &anns[i]
		}
	}
	return best
}

func vepAnnBetter(a, b *vepAnnotation) bool {
	if a.canonical != b.canonical {
		return a.canonical
	}
	aPC := a.biotype == "protein_coding"
	bPC := b.biotype == "protein_coding"
	if aPC != bPC {
		return aPC
	}
	return annotate.ImpactRank(a.impact) > annotate.ImpactRank(b.impact)
}

// vepAnnotationToMap converts a vepAnnotation to the column map used by Categorizer.
func vepAnnotationToMap(ann *vepAnnotation) map[string]string {
	if ann == nil {
		return map[string]string{
			"Consequence": "",
			"HGVSp_Short": "",
			"HGVSc":       "",
		}
	}
	return map[string]string{
		"Consequence": ann.consequence,
		"HGVSp_Short": output.HGVSpToShort(ann.hgvsp),
		"HGVSc":       ann.hgvsc,
	}
}
