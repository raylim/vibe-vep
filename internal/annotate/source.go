package annotate

import "github.com/inodb/vibe-vep/internal/vcf"

// MatchLevel describes what a source matches on.
type MatchLevel string

const (
	// MatchGenomic matches on chr:pos:ref:alt (genomic coordinates).
	// Sources using this level are assembly-specific (GRCh37 vs GRCh38).
	MatchGenomic MatchLevel = "genomic"
	// MatchProteinPosition matches on transcript + amino acid position.
	// Sources using this level are transcript-version sensitive.
	MatchProteinPosition MatchLevel = "protein_position"
	// MatchGene matches on gene symbol only.
	MatchGene MatchLevel = "gene"
)

// AnnotationSource adds external data to variant annotations.
type AnnotationSource interface {
	Name() string         // e.g. "alphamissense"
	Version() string      // e.g. "2023"
	MatchLevel() MatchLevel // what the source matches on
	Columns() []ColumnDef // columns this source provides
	Annotate(v *vcf.Variant, anns []*Annotation)
}

// ColumnDef describes a column provided by an annotation source.
type ColumnDef struct {
	Name        string // short name, e.g. "score"
	Description string // human-readable description
}

// CoreColumns defines the columns produced by vibe-vep's core prediction.
var CoreColumns = []ColumnDef{
	{Name: "hugo_symbol", Description: "Gene symbol"},
	{Name: "consequence", Description: "SO consequence term"},
	{Name: "variant_classification", Description: "MAF variant classification"},
	{Name: "transcript_id", Description: "Ensembl transcript ID"},
	{Name: "hgvsc", Description: "HGVS coding DNA notation"},
	{Name: "hgvsp", Description: "HGVS protein notation (3-letter)"},
	{Name: "hgvsp_short", Description: "HGVS protein notation (1-letter)"},
	{Name: "canonical_mskcc", Description: "MSK canonical transcript"},
	{Name: "canonical_ensembl", Description: "Ensembl canonical transcript"},
}

// SetExtra sets a value in the annotation's Extra map.
func (a *Annotation) SetExtra(source, field, value string) {
	a.SetExtraKey(source+"."+field, value)
}

// SetExtraKey sets a value in the annotation's Extra map using a pre-built key.
func (a *Annotation) SetExtraKey(key, value string) {
	if a.Extra == nil {
		a.Extra = make(map[string]string)
	}
	a.Extra[key] = value
}

// GetExtra returns a value from the annotation's Extra map.
func (a *Annotation) GetExtra(source, field string) string {
	return a.GetExtraKey(source + "." + field)
}

// GetExtraKey returns a value from the annotation's Extra map using a pre-built key.
func (a *Annotation) GetExtraKey(key string) string {
	if a.Extra == nil {
		return ""
	}
	return a.Extra[key]
}
