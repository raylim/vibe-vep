package hotspots

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const testTSV = `hugo_symbol	residue	reference_amino_acid	amino_acid_position	cluster	pdb_chains	class	variant_amino_acid	q_value	qvalue_pancan	qvaluect	p_value	tumor_count	tumor_type_count	tumor_type_composition	indel_size	transcript_id	type	missense_count	trunc_count	inframe_count	splice_count	total_count	missense_fraction	trunc_fraction	inframe_fraction	splice_fraction
BRAF	V600	V:897	600				R:4|E:833	0.0	0.0	0.0		897		skin:357		ENST00000288602	single residue	897.0	0.0	0.0	0.0	897.0	1.0	0.0	0.0	0.0
KRAS	G12	G:2175	12				A:134	0.0	0.0	0.0		2175		lung:570		ENST00000256078	single residue	2175.0	0.0	0.0	0.0	2175.0	1.0	0.0	0.0	0.0
KRAS	G13	G:264	13				D:214	0.0	0.0	5.49e-287		264		bowel:130		ENST00000256078	single residue	264.0	0.0	0.0	0.0	264.0	1.0	0.0	0.0	0.0
TP53	R175	R:416	175				H:386	0.0	0.0	1.74e-111		416		bowel:108		ENST00000269305	single residue	416.0	0.0	0.0	0.0	416.0	1.0	0.0	0.0	0.0
EGFR	L858	L:144	858				R:144	0.0	3.41e-276	0.0		144		lung:143		ENST00000275493	single residue	144.0	0.0	0.0	0.0	144.0	1.0	0.0	0.0	0.0
`

func writeTSV(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "hotspots.txt")
	require.NoError(t, os.WriteFile(path, []byte(testTSV), 0644))
	return path
}

func TestLoadAndLookup(t *testing.T) {
	store, err := Load(writeTSV(t))
	require.NoError(t, err)

	assert.Equal(t, 4, store.GeneCount()) // BRAF, KRAS, TP53, EGFR
	assert.Equal(t, 5, store.HotspotCount())

	// BRAF V600 is a hotspot
	h, ok := store.Lookup("BRAF", 600)
	assert.True(t, ok)
	assert.Equal(t, "single residue", h.Type)
	assert.Equal(t, 0.0, h.QValue)

	// KRAS G12 is a hotspot
	h, ok = store.Lookup("KRAS", 12)
	assert.True(t, ok)
	assert.Equal(t, int64(12), h.Position)

	// KRAS G13 is a hotspot
	_, ok = store.Lookup("KRAS", 13)
	assert.True(t, ok)

	// KRAS position 14 is NOT a hotspot
	_, ok = store.Lookup("KRAS", 14)
	assert.False(t, ok)

	// Unknown gene
	_, ok = store.Lookup("FAKEGENE", 100)
	assert.False(t, ok)
}

func TestAnnotate(t *testing.T) {
	store, err := Load(writeTSV(t))
	require.NoError(t, err)
	src := NewSource(store)

	assert.Equal(t, "hotspots", src.Name())
	assert.Equal(t, "v2", src.Version())
	assert.Len(t, src.Columns(), 3)

	v := &vcf.Variant{Chrom: "7", Pos: 140753336, Ref: "A", Alt: "T"}
	anns := []*annotate.Annotation{
		{GeneName: "BRAF", ProteinPosition: 600, Consequence: "missense_variant"},
		{GeneName: "BRAF", ProteinPosition: 100, Consequence: "missense_variant"}, // not a hotspot
		{GeneName: "", ProteinPosition: 0, Consequence: "intergenic_variant"},      // no gene
	}

	src.Annotate(v, anns)

	assert.Equal(t, "Y", anns[0].GetExtraKey(extraKeyHotspot))
	assert.Equal(t, "single residue", anns[0].GetExtraKey(extraKeyType))
	assert.Equal(t, "0", anns[0].GetExtraKey(extraKeyQValue))

	assert.Equal(t, "", anns[1].GetExtraKey(extraKeyHotspot)) // not a hotspot
	assert.Equal(t, "", anns[2].GetExtraKey(extraKeyHotspot)) // no gene
}

func TestFormatQValue(t *testing.T) {
	assert.Equal(t, "0", FormatQValue(0))
	assert.Equal(t, "0.05", FormatQValue(0.05))
	assert.Equal(t, "5.490e-287", FormatQValue(5.49e-287))
}
