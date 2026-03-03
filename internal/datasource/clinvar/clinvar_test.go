package clinvar

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const testVCF = `##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease name">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
7	140753336	846933	T	A	.	.	ALLELEID=826271;CLNDN=Melanoma;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel
12	25245350	RCV000038388	C	A	.	.	ALLELEID=37050;CLNDN=Lung_adenocarcinoma;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter
17	7675088	RCV000013130	C	T	.	.	ALLELEID=22222;CLNDN=Li-Fraumeni_syndrome;CLNSIG=Pathogenic/Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel
1	100	.	A	G	.	.	ALLELEID=99;CLNDN=not_provided;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts
`

func writeVCF(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "clinvar.vcf")
	require.NoError(t, os.WriteFile(path, []byte(testVCF), 0644))
	return path
}

func TestLoadAndLookup(t *testing.T) {
	store, err := Load(writeVCF(t))
	require.NoError(t, err)
	assert.Equal(t, 4, store.Count())

	// BRAF V600E — Pathogenic
	e, ok := store.Lookup("7", 140753336, "T", "A")
	assert.True(t, ok)
	assert.Equal(t, "Pathogenic", e.ClnSig)
	assert.Equal(t, "reviewed_by_expert_panel", e.RevStat)
	assert.Equal(t, "Melanoma", e.ClnDN)

	// KRAS — Pathogenic
	e, ok = store.Lookup("12", 25245350, "C", "A")
	assert.True(t, ok)
	assert.Equal(t, "Pathogenic", e.ClnSig)

	// With chr prefix
	e, ok = store.Lookup("chr7", 140753336, "T", "A")
	assert.True(t, ok)
	assert.Equal(t, "Pathogenic", e.ClnSig)

	// Wrong alt
	_, ok = store.Lookup("7", 140753336, "T", "C")
	assert.False(t, ok)

	// Unknown position
	_, ok = store.Lookup("7", 999999999, "A", "T")
	assert.False(t, ok)
}

func TestGobRoundTrip(t *testing.T) {
	store, err := Load(writeVCF(t))
	require.NoError(t, err)

	gobPath := filepath.Join(t.TempDir(), "clinvar.gob")
	require.NoError(t, store.SaveGob(gobPath))

	store2, err := LoadGob(gobPath)
	require.NoError(t, err)
	assert.Equal(t, store.Count(), store2.Count())

	// Verify data integrity after round-trip
	e, ok := store2.Lookup("7", 140753336, "T", "A")
	assert.True(t, ok)
	assert.Equal(t, "Pathogenic", e.ClnSig)
	assert.Equal(t, "Melanoma", e.ClnDN)
}

func TestAnnotate(t *testing.T) {
	store, err := Load(writeVCF(t))
	require.NoError(t, err)
	src := NewSource(store)

	assert.Equal(t, "clinvar", src.Name())
	assert.Len(t, src.Columns(), 3)

	v := &vcf.Variant{Chrom: "7", Pos: 140753336, Ref: "T", Alt: "A"}
	anns := []*annotate.Annotation{
		{GeneName: "BRAF", Consequence: "missense_variant"},
		{GeneName: "BRAF", Consequence: "downstream_gene_variant"},
	}

	src.Annotate(v, anns)

	// Both annotations get ClinVar data (variant-level)
	for _, ann := range anns {
		assert.Equal(t, "Pathogenic", ann.GetExtraKey(extraKeyClnSig))
		assert.Equal(t, "reviewed_by_expert_panel", ann.GetExtraKey(extraKeyRevStat))
		assert.Equal(t, "Melanoma", ann.GetExtraKey(extraKeyClnDN))
	}

	// Variant not in ClinVar
	v2 := &vcf.Variant{Chrom: "7", Pos: 999999, Ref: "A", Alt: "T"}
	anns2 := []*annotate.Annotation{{GeneName: "X"}}
	src.Annotate(v2, anns2)
	assert.Equal(t, "", anns2[0].GetExtraKey(extraKeyClnSig))
}

func TestExtractInfo(t *testing.T) {
	info := "ALLELEID=826271;CLNDN=Melanoma;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel"
	assert.Equal(t, "Pathogenic", extractInfo(info, "CLNSIG="))
	assert.Equal(t, "Melanoma", extractInfo(info, "CLNDN="))
	assert.Equal(t, "reviewed_by_expert_panel", extractInfo(info, "CLNREVSTAT="))
	assert.Equal(t, "", extractInfo(info, "MISSING="))
}
