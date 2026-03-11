package output_test

// ANNOVAR annotation loader for the ClinVar benchmark.
//
// ANNOVAR outputs a tab-delimited .hg38_multianno.txt file with one row per
// variant.  The AAChange column contains comma-separated transcript entries of
// the form GENE:NM_ID:exonN:c.HGVSc:p.HGVSp, where protein HGVS uses
// single-letter amino acid codes (e.g. p.T141A).
//
// This file provides:
//   - convertANNOVARProtein  — 1-letter → 3-letter AA conversion
//   - annovarVariantMap      — lookup by exact key or genomic position
//   - loadAnnovarTxt         — parse the multianno.txt file
//   - annovarConsequenceToSO — map ANNOVAR exonic function to SO terms
//   - pickBestAnnovarByImpact

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
	"unicode"

	"github.com/inodb/vibe-vep/internal/vcf"
)

// aa1to3 maps single-letter amino acid codes to three-letter codes.
var aa1to3 = map[byte]string{
	'A': "Ala", 'R': "Arg", 'N': "Asn", 'D': "Asp", 'C': "Cys",
	'E': "Glu", 'Q': "Gln", 'G': "Gly", 'H': "His", 'I': "Ile",
	'L': "Leu", 'K': "Lys", 'M': "Met", 'F': "Phe", 'P': "Pro",
	'S': "Ser", 'T': "Thr", 'W': "Trp", 'Y': "Tyr", 'V': "Val",
	'X': "Ter", // ANNOVAR uses X for stop codons
}

// isAA1 returns true if ch is a valid single-letter amino acid code (uppercase).
func isAA1(ch byte) bool {
	_, ok := aa1to3[ch]
	return ok
}

// convertANNOVARProtein converts ANNOVAR 1-letter amino acid codes in a
// protein HGVS string to 3-letter codes.
//
// ANNOVAR examples:
//   - p.T141A        → p.Thr141Ala
//   - p.R793X        → p.Arg793Ter
//   - p.L298Rfs*17   → p.Leu298Argfs*17
//   - p.K237_I238delinsN → p.Lys237_Ile238delinsAsn
//
// The algorithm walks the string and converts uppercase amino acid letters
// when they appear adjacent to digits or in positions that indicate they
// are amino acid codes rather than part of a keyword like "del", "ins", "fs".
func convertANNOVARProtein(p string) string {
	if p == "" || p == "." || p == "UNKNOWN" {
		return p
	}
	prefix := ""
	body := p
	if strings.HasPrefix(p, "p.") {
		prefix = "p."
		body = p[2:]
	}

	var out strings.Builder
	out.WriteString(prefix)

	n := len(body)
	for i := 0; i < n; {
		ch := body[i]

		// Only consider uppercase letters for AA conversion.
		if ch >= 'A' && ch <= 'Z' && isAA1(ch) {
			prevIsDigit := i > 0 && unicode.IsDigit(rune(body[i-1]))
			nextIsDigit := i+1 < n && unicode.IsDigit(rune(body[i+1]))
			atEnd := i+1 == n
			nextIsUpperAA := i+1 < n && body[i+1] >= 'A' && body[i+1] <= 'Z' && isAA1(body[i+1])
			prevIsUnderscore := i > 0 && body[i-1] == '_'
			// Don't convert if previous char is lowercase — mid-word in a keyword
			// like "ins", "del", "dup", "fs".
			prevIsLower := i > 0 && body[i-1] >= 'a' && body[i-1] <= 'z'

			shouldConvert := nextIsDigit || prevIsDigit || atEnd || nextIsUpperAA || prevIsUnderscore
			if prevIsLower {
				shouldConvert = false
			}

			if shouldConvert {
				out.WriteString(aa1to3[ch])
				i++
				continue
			}
		}

		out.WriteByte(ch)
		i++
	}

	return out.String()
}

// annovarAnnotation holds one transcript annotation from ANNOVAR's AAChange column.
type annovarAnnotation struct {
	geneSymbol string
	transcript string // NM_ accession (may include version)
	hgvsc      string
	hgvsp      string // converted to 3-letter AA codes
	exonicFunc string // from ExonicFunc column
}

// annovarVariantMap maps variant positions to ANNOVAR annotations.
type annovarVariantMap struct {
	// exact: "chr:pos:ref:alt" → annotations
	exact map[string][]annovarAnnotation
	// byPos: "chr:pos" → annotations (fallback for indels with shifted coordinates)
	byPos map[string][]annovarAnnotation
}

func newAnnovarVariantMap() annovarVariantMap {
	return annovarVariantMap{
		exact: make(map[string][]annovarAnnotation),
		byPos: make(map[string][]annovarAnnotation),
	}
}

// register stores annotations under exact and positional keys.
// For ANNOVAR deletions (Alt == "-"), Start = VCF POS+1 (anchor stripped),
// so we also register under pos-1 to match VCF-based positional lookups.
func (m *annovarVariantMap) register(chrom string, pos int, ref, alt string, anns []annovarAnnotation) {
	isAnnovarDel := alt == "-"

	exactKey := snpEffExactKey(normalizeChrom(chrom), int64(pos), ref, alt)
	m.exact[exactKey] = append(m.exact[exactKey], anns...)

	posKey := snpEffPosKey(normalizeChrom(chrom), int64(pos))
	m.byPos[posKey] = append(m.byPos[posKey], anns...)

	if isAnnovarDel {
		// Also register under VCF POS (pos-1) for positional lookup.
		posKeyMinus1 := snpEffPosKey(normalizeChrom(chrom), int64(pos-1))
		m.byPos[posKeyMinus1] = append(m.byPos[posKeyMinus1], anns...)
	}
}

// lookup returns ANNOVAR annotations for a vcf.Variant, trying exact match first.
func (m annovarVariantMap) lookup(v *vcf.Variant) []annovarAnnotation {
	key := snpEffExactKey(v.NormalizeChrom(), v.Pos, v.Ref, v.Alt)
	if anns, ok := m.exact[key]; ok {
		return anns
	}
	return m.byPos[snpEffPosKey(v.NormalizeChrom(), v.Pos)]
}

// loadAnnovarTxt parses an ANNOVAR .hg38_multianno.txt file (plain or gzipped).
func loadAnnovarTxt(path string) (annovarVariantMap, error) {
	f, err := os.Open(path)
	if err != nil {
		return annovarVariantMap{}, fmt.Errorf("open %s: %w", path, err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gr, err := gzip.NewReader(f)
		if err != nil {
			return annovarVariantMap{}, fmt.Errorf("gzip reader: %w", err)
		}
		defer gr.Close()
		r = gr
	}

	m := newAnnovarVariantMap()
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 4*1024*1024), 4*1024*1024)

	// Parse header to locate column indices.
	if !sc.Scan() {
		return m, nil
	}
	header := strings.Split(sc.Text(), "\t")

	idxOf := func(name string) int {
		for i, h := range header {
			if strings.EqualFold(h, name) {
				return i
			}
		}
		return -1
	}

	iChr := idxOf("Chr")
	iStart := idxOf("Start")
	iRef := idxOf("Ref")
	iAlt := idxOf("Alt")

	iExonicFunc := -1
	iAAChange := -1
	for i, h := range header {
		lh := strings.ToLower(h)
		if strings.HasPrefix(lh, "exonicfunc.") && iExonicFunc < 0 {
			iExonicFunc = i
		}
		if strings.HasPrefix(lh, "aachange.") && iAAChange < 0 {
			iAAChange = i
		}
	}

	if iChr < 0 || iStart < 0 || iRef < 0 || iAlt < 0 {
		return m, fmt.Errorf("missing required columns in ANNOVAR header: %v", header)
	}

	safeField := func(fields []string, i int) string {
		if i < 0 || i >= len(fields) {
			return "."
		}
		return fields[i]
	}

	for sc.Scan() {
		fields := strings.Split(sc.Text(), "\t")
		if len(fields) <= iChr {
			continue
		}

		chrom := fields[iChr]
		var pos int
		fmt.Sscanf(safeField(fields, iStart), "%d", &pos)
		ref := safeField(fields, iRef)
		alt := safeField(fields, iAlt)
		exonicFunc := safeField(fields, iExonicFunc)
		aaChange := safeField(fields, iAAChange)

		if aaChange == "." || aaChange == "" {
			continue
		}

		anns := parseAAChange(aaChange, exonicFunc)
		if len(anns) == 0 {
			continue
		}

		m.register(chrom, pos, ref, alt, anns)
	}
	if err := sc.Err(); err != nil {
		return m, fmt.Errorf("scan %s: %w", path, err)
	}
	return m, nil
}

// parseAAChange parses the AAChange.refGeneWithVer column, which contains
// comma-separated entries like "GENE:NM_ID:exonN:c.xxx:p.xxx".
func parseAAChange(aaChange, exonicFunc string) []annovarAnnotation {
	var anns []annovarAnnotation
	for _, entry := range strings.Split(aaChange, ",") {
		entry = strings.TrimSpace(entry)
		if entry == "" || entry == "." || entry == "UNKNOWN" {
			continue
		}
		parts := strings.SplitN(entry, ":", 5)
		if len(parts) < 2 {
			continue
		}
		ann := annovarAnnotation{
			geneSymbol: parts[0],
			transcript: parts[1],
			exonicFunc: exonicFunc,
		}
		if len(parts) >= 4 {
			ann.hgvsc = parts[3]
		}
		if len(parts) >= 5 {
			ann.hgvsp = convertANNOVARProtein(parts[4])
		}
		anns = append(anns, ann)
	}
	return anns
}

// annovarConsequenceToSO maps ANNOVAR's ExonicFunc values to SO consequence terms.
func annovarConsequenceToSO(exonicFunc string) string {
	switch strings.ToLower(strings.TrimSpace(exonicFunc)) {
	case "nonsynonymous snv":
		return "missense_variant"
	case "stopgain":
		return "stop_gained"
	case "stoploss":
		return "stop_lost"
	case "startloss":
		return "start_lost"
	case "synonymous snv":
		return "synonymous_variant"
	case "frameshift deletion", "frameshift insertion", "frameshift substitution":
		return "frameshift_variant"
	case "nonframeshift deletion":
		return "inframe_deletion"
	case "nonframeshift insertion":
		return "inframe_insertion"
	case "nonframeshift substitution":
		return "inframe_insertion"
	}
	return exonicFunc
}

// annovarImpactRank returns a numeric impact rank for an ANNOVAR exonic function.
func annovarImpactRank(exonicFunc string) int {
	switch strings.ToLower(strings.TrimSpace(exonicFunc)) {
	case "stopgain", "stoploss", "startloss",
		"frameshift deletion", "frameshift insertion", "frameshift substitution":
		return 4
	case "nonsynonymous snv":
		return 3
	case "nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution":
		return 2
	case "synonymous snv":
		return 1
	}
	return 0
}

// pickBestAnnovarByImpact returns the highest-impact ANNOVAR annotation
// that has a non-empty HGVSp.
func pickBestAnnovarByImpact(anns []annovarAnnotation) *annovarAnnotation {
	var best *annovarAnnotation
	for i := range anns {
		a := &anns[i]
		if a.hgvsp == "" {
			continue
		}
		if best == nil || annovarImpactRank(a.exonicFunc) > annovarImpactRank(best.exonicFunc) {
			best = a
		}
	}
	return best
}
