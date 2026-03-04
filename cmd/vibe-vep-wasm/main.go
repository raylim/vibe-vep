//go:build js && wasm

// Command vibe-vep-wasm provides a WebAssembly build of the vibe-vep annotation
// engine for use in browser-based tutorials.
package main

import (
	"embed"
	"encoding/json"
	"fmt"
	"sort"
	"strings"
	"syscall/js"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/maf"
	"github.com/inodb/vibe-vep/internal/vcf"
)

//go:embed data/transcripts_12.json data/transcripts_7.json
var transcriptFS embed.FS

//go:embed data/example.vcf
var exampleVCF string

//go:embed data/example.maf
var exampleMAF string

var (
	txCache   *cache.Cache
	ann       *annotate.Annotator
	wasmReady bool
)

func main() {
	if err := initialize(); err != nil {
		js.Global().Get("console").Call("error", "vibe-vep init failed: "+err.Error())
		return
	}

	js.Global().Set("vibeVEP", map[string]interface{}{
		"annotateVariant": js.FuncOf(annotateVariant),
		"annotateVCF":     js.FuncOf(annotateVCF),
		"annotateMAF":     js.FuncOf(annotateMAF),
		"help":            js.FuncOf(help),
		"version":         js.FuncOf(version),
		"listGenes":       js.FuncOf(listGenes),
		"exampleVCF":      js.FuncOf(getExampleVCF),
		"exampleMAF":      js.FuncOf(getExampleMAF),
		"ready":           js.FuncOf(ready),
	})

	// Signal that WASM is ready.
	js.Global().Get("console").Call("log", "vibe-vep WASM ready")
	if cb := js.Global().Get("onVibeVEPReady"); !cb.IsUndefined() {
		cb.Invoke()
	}

	// Block forever so the Go runtime stays alive.
	select {}
}

func initialize() error {
	txCache = cache.New()

	files := []string{"data/transcripts_12.json", "data/transcripts_7.json"}
	for _, f := range files {
		data, err := transcriptFS.ReadFile(f)
		if err != nil {
			return fmt.Errorf("read %s: %w", f, err)
		}
		var transcripts []*cache.Transcript
		if err := json.Unmarshal(data, &transcripts); err != nil {
			return fmt.Errorf("decode %s: %w", f, err)
		}
		for _, t := range transcripts {
			txCache.AddTranscript(t)
		}
	}

	txCache.BuildIndex()
	ann = annotate.NewAnnotator(txCache)
	wasmReady = true
	return nil
}

func annotateVariant(_ js.Value, args []js.Value) interface{} {
	defer func() {
		if r := recover(); r != nil {
			js.Global().Get("console").Call("error", fmt.Sprintf("panic: %v", r))
		}
	}()

	if len(args) < 1 {
		return "Error: missing variant specification"
	}
	specInput := args[0].String()

	spec, err := annotate.ParseVariantSpec(specInput)
	if err != nil {
		return "Error: " + err.Error()
	}

	var variants []*vcf.Variant
	var header strings.Builder

	switch spec.Type {
	case annotate.SpecGenomic:
		variants = []*vcf.Variant{{
			Chrom: spec.Chrom,
			Pos:   spec.Pos,
			Ref:   spec.Ref,
			Alt:   spec.Alt,
		}}

	case annotate.SpecProtein:
		variants, err = annotate.ReverseMapProteinChange(txCache, spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
		if err != nil {
			return "Error: " + err.Error()
		}
		fmt.Fprintf(&header, "Query: %s %c%d%c\n", spec.GeneName, spec.RefAA, spec.Position, spec.AltAA)
		fmt.Fprintf(&header, "Found %d genomic variant(s):\n", len(variants))
		for _, v := range variants {
			fmt.Fprintf(&header, "  %s:%d %s>%s\n", v.Chrom, v.Pos, v.Ref, v.Alt)
		}
		header.WriteString("\n")

	case annotate.SpecHGVSc:
		variants, err = annotate.ReverseMapHGVSc(txCache, spec.TranscriptID, spec.CDSChange)
		if err != nil {
			return "Error: " + err.Error()
		}
		fmt.Fprintf(&header, "Query: %s c.%s\n", spec.TranscriptID, spec.CDSChange)
		fmt.Fprintf(&header, "Found %d genomic variant(s):\n", len(variants))
		for _, v := range variants {
			fmt.Fprintf(&header, "  %s:%d %s>%s\n", v.Chrom, v.Pos, v.Ref, v.Alt)
		}
		header.WriteString("\n")
	}

	var out strings.Builder
	out.WriteString(header.String())

	for _, v := range variants {
		fmt.Fprintf(&out, "Variant: %s:%d %s>%s\n\n", v.Chrom, v.Pos, v.Ref, v.Alt)

		anns, err := ann.Annotate(v)
		if err != nil {
			fmt.Fprintf(&out, "Error: %s\n", err)
			continue
		}

		out.WriteString(formatAnnotationTable(anns))
	}

	return out.String()
}

func annotateVCF(_ js.Value, args []js.Value) interface{} {
	defer func() {
		if r := recover(); r != nil {
			js.Global().Get("console").Call("error", fmt.Sprintf("panic: %v", r))
		}
	}()

	if len(args) < 1 {
		return "Error: missing VCF content"
	}
	content := args[0].String()

	parser, err := vcf.NewParserFromReader(strings.NewReader(content))
	if err != nil {
		return "Error: " + err.Error()
	}

	var out strings.Builder
	count := 0

	for {
		v, err := parser.Next()
		if err != nil {
			fmt.Fprintf(&out, "Error at variant %d: %s\n", count+1, err)
			break
		}
		if v == nil {
			break
		}
		count++

		anns, err := ann.Annotate(v)
		if err != nil {
			fmt.Fprintf(&out, "Error annotating %s:%d: %s\n", v.Chrom, v.Pos, err)
			continue
		}

		fmt.Fprintf(&out, "Variant %d: %s:%d %s>%s\n\n", count, v.Chrom, v.Pos, v.Ref, v.Alt)
		out.WriteString(formatAnnotationTable(anns))
		out.WriteString("\n")
	}

	if count == 0 {
		return "No variants found in VCF content."
	}

	fmt.Fprintf(&out, "Annotated %d variant(s).\n", count)
	return out.String()
}

func annotateMAF(_ js.Value, args []js.Value) interface{} {
	defer func() {
		if r := recover(); r != nil {
			js.Global().Get("console").Call("error", fmt.Sprintf("panic: %v", r))
		}
	}()

	if len(args) < 1 {
		return "Error: missing MAF content"
	}
	content := args[0].String()

	parser, err := maf.NewParserFromReader(strings.NewReader(content))
	if err != nil {
		return "Error: " + err.Error()
	}

	var out strings.Builder
	count := 0

	for {
		v, err := parser.Next()
		if err != nil {
			fmt.Fprintf(&out, "Error at variant %d: %s\n", count+1, err)
			break
		}
		if v == nil {
			break
		}
		count++

		anns, err := ann.Annotate(v)
		if err != nil {
			fmt.Fprintf(&out, "Error annotating %s:%d: %s\n", v.Chrom, v.Pos, err)
			continue
		}

		fmt.Fprintf(&out, "Variant %d: %s:%d %s>%s\n\n", count, v.Chrom, v.Pos, v.Ref, v.Alt)
		out.WriteString(formatAnnotationTable(anns))
		out.WriteString("\n")
	}

	if count == 0 {
		return "No variants found in MAF content."
	}

	fmt.Fprintf(&out, "Annotated %d variant(s).\n", count)
	return out.String()
}

func help(_ js.Value, _ []js.Value) interface{} {
	return `vibe-vep — variant effect predictor (WASM demo)

Commands:
  vibe-vep annotate variant <spec>   Annotate a single variant
  vibe-vep annotate vcf <file>       Annotate variants from a VCF file
  vibe-vep annotate maf <file>       Annotate variants from a MAF file
  vibe-vep version                   Show version info
  vibe-vep genes                     List available genes

Variant specification formats:
  Genomic:  12:25245351:C:A
  Protein:  KRAS G12C
  HGVSc:    KRAS c.35G>T

Shell commands:
  cat <file>    Show file contents
  ls            List available files
  clear         Clear terminal
  help          Show this help

This demo includes KRAS and EGFR transcripts (GRCh38).
`
}

func version(_ js.Value, _ []js.Value) interface{} {
	chroms := txCache.Chromosomes()
	genes := listAllGenes()
	return fmt.Sprintf(`vibe-vep (WASM demo)
Assembly: GRCh38
Transcripts: %d
Chromosomes: %s
Genes: %s
`, txCache.TranscriptCount(), strings.Join(chroms, ", "), strings.Join(genes, ", "))
}

func listGenes(_ js.Value, _ []js.Value) interface{} {
	genes := listAllGenes()
	return strings.Join(genes, "\n")
}

func getExampleVCF(_ js.Value, _ []js.Value) interface{} {
	return exampleVCF
}

func getExampleMAF(_ js.Value, _ []js.Value) interface{} {
	return exampleMAF
}

func ready(_ js.Value, _ []js.Value) interface{} {
	return wasmReady
}

func listAllGenes() []string {
	seen := make(map[string]bool)
	for _, chrom := range txCache.Chromosomes() {
		for _, t := range txCache.FindTranscriptsByChrom(chrom) {
			if t.GeneName != "" {
				seen[t.GeneName] = true
			}
		}
	}
	genes := make([]string, 0, len(seen))
	for g := range seen {
		genes = append(genes, g)
	}
	sort.Strings(genes)
	return genes
}

func formatAnnotationTable(anns []*annotate.Annotation) string {
	var out strings.Builder

	// Calculate column widths for alignment.
	type row struct {
		gene, tx, canon, csq, impact, hgvsc, hgvsp string
	}
	rows := make([]row, len(anns))
	widths := [7]int{4, 10, 5, 11, 6, 5, 5} // header widths

	for i, a := range anns {
		canon := "no"
		if a.IsCanonical {
			canon = "YES"
		}
		r := row{a.GeneName, a.TranscriptID, canon, a.Consequence, a.Impact, a.HGVSc, a.HGVSp}
		rows[i] = r
		for j, s := range []string{r.gene, r.tx, r.canon, r.csq, r.impact, r.hgvsc, r.hgvsp} {
			if len(s) > widths[j] {
				widths[j] = len(s)
			}
		}
	}

	fmtStr := fmt.Sprintf("%%-%ds  %%-%ds  %%-%ds  %%-%ds  %%-%ds  %%-%ds  %%s\n",
		widths[0], widths[1], widths[2], widths[3], widths[4], widths[5])

	fmt.Fprintf(&out, fmtStr, "Gene", "Transcript", "Canon", "Consequence", "Impact", "HGVSc", "HGVSp")
	for _, r := range rows {
		fmt.Fprintf(&out, fmtStr, r.gene, r.tx, r.canon, r.csq, r.impact, r.hgvsc, r.hgvsp)
	}

	return out.String()
}
