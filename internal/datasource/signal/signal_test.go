package signal

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const testTSV = `Hugo_Symbol	Chromosome	Start_Position	End_Position	Variant_Classification	Variant_Type	Strand	Reference_Allele	Alternate_Allele	n_impact	f_biallelic	n_germline_homozygous	n_eur	n_afr	n_asn	n_asj	n_oth	f_impact	f_eur	f_afr	f_asn	f_asj	f_oth
BRCA2	13	32890572	32890572	Missense_Mutation	SNP	+	G	A	5	0	0	3	1	0	1	0	0.000291511	0.000305064	0.00084034	0	0.000318877	0
BRCA1	17	43093220	43093220	Nonsense_Mutation	SNP	+	C	T	3	1	0	0	0	0	3	0	0.000174907	0	0	0	0.000956633	0
ATM	11	108288900	108288900	Missense_Mutation	SNP	+	C	T	12	0.5	0	8	1	0	2	1	0.000699627	0.000813504	0.00084034	0	0.000637755	0.000475285
APC	5	112175770	112175771	Frame_Shift_Del	DEL	+	AG	-	2	1	0	2	0	0	0	0	0.000116604	0.000203376	0	0	0	0
`

func writeTSV(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "signal.txt")
	require.NoError(t, os.WriteFile(path, []byte(testTSV), 0644))
	return path
}

func TestLoadAndLookup(t *testing.T) {
	store, err := Load(writeTSV(t))
	require.NoError(t, err)
	assert.Equal(t, 4, store.Count())

	// BRCA2 variant
	e, ok := store.Lookup("13", 32890572, "G", "A")
	assert.True(t, ok)
	assert.Equal(t, "BRCA2", e.Gene)
	assert.Equal(t, 5, e.CountAll)
	assert.InDelta(t, 0.000291511, e.FreqAll, 1e-9)

	// APC deletion (- converted to empty)
	e, ok = store.Lookup("5", 112175770, "AG", "")
	assert.True(t, ok)
	assert.Equal(t, "APC", e.Gene)

	// With chr prefix
	e, ok = store.Lookup("chr13", 32890572, "G", "A")
	assert.True(t, ok)
	assert.Equal(t, "BRCA2", e.Gene)

	// Unknown variant
	_, ok = store.Lookup("13", 999999999, "A", "T")
	assert.False(t, ok)
}

func TestGobRoundTrip(t *testing.T) {
	store, err := Load(writeTSV(t))
	require.NoError(t, err)

	gobPath := filepath.Join(t.TempDir(), "signal.gob")
	require.NoError(t, store.SaveGob(gobPath))

	store2, err := LoadGob(gobPath)
	require.NoError(t, err)
	assert.Equal(t, store.Count(), store2.Count())

	e, ok := store2.Lookup("13", 32890572, "G", "A")
	assert.True(t, ok)
	assert.Equal(t, "BRCA2", e.Gene)
	assert.Equal(t, 5, e.CountAll)
}

func TestAnnotate(t *testing.T) {
	store, err := Load(writeTSV(t))
	require.NoError(t, err)
	src := NewSource(store)

	assert.Equal(t, "signal", src.Name())
	assert.Len(t, src.Columns(), 3)

	v := &vcf.Variant{Chrom: "13", Pos: 32890572, Ref: "G", Alt: "A"}
	anns := []*annotate.Annotation{
		{GeneName: "BRCA2", Consequence: "missense_variant"},
	}

	src.Annotate(v, anns)
	assert.Equal(t, "germline", anns[0].GetExtraKey(extraKeyMutationStatus))
	assert.Equal(t, "5", anns[0].GetExtraKey(extraKeyCountCarriers))
	assert.NotEmpty(t, anns[0].GetExtraKey(extraKeyFrequency))

	// Variant not in SIGNAL
	v2 := &vcf.Variant{Chrom: "1", Pos: 100, Ref: "A", Alt: "T"}
	anns2 := []*annotate.Annotation{{GeneName: "X"}}
	src.Annotate(v2, anns2)
	assert.Equal(t, "", anns2[0].GetExtraKey(extraKeyMutationStatus))
}
