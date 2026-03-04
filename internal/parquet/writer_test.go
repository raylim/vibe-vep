package parquet

import (
	"os"
	"testing"

	"github.com/inodb/vibe-vep/internal/annotate"
	pq "github.com/parquet-go/parquet-go"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestChromToNumeric(t *testing.T) {
	tests := []struct {
		chrom string
		want  int32
	}{
		{"1", 1},
		{"12", 12},
		{"22", 22},
		{"X", 23},
		{"Y", 24},
		{"MT", 25},
		{"chr1", 1},
		{"chrX", 23},
		{"chrMT", 25},
		{"chrM", 25},
		{"unknown", 0},
	}
	for _, tt := range tests {
		t.Run(tt.chrom, func(t *testing.T) {
			assert.Equal(t, tt.want, ChromToNumeric(tt.chrom))
		})
	}
}

func TestRoundTrip(t *testing.T) {
	rows := []Row{
		{
			ChromNumeric: 12, Pos: 25245350, Chrom: "12",
			Ref: "C", Alt: "A",
			TranscriptID: "ENST00000311936.8", GeneName: "KRAS", GeneID: "ENSG00000133703",
			Consequence: "missense_variant", Impact: "MODERATE",
			CDSPosition: 35, ProteinPosition: 12, AminoAcidChange: "G12C",
			IsCanonical: true, Biotype: "protein_coding",
			HGVSp: "p.Gly12Cys", HGVSc: "c.34G>T",
			AMScore: 0.9876, AMClass: "likely_pathogenic",
			OncokbGeneType: "ONCOGENE",
		},
		{
			ChromNumeric: 1, Pos: 115256530, Chrom: "1",
			Ref: "T", Alt: "C",
			TranscriptID: "ENST00000369535.9", GeneName: "NRAS", GeneID: "ENSG00000213281",
			Consequence: "synonymous_variant", Impact: "LOW",
			CDSPosition: 60, ProteinPosition: 20,
			IsCanonical: true, Biotype: "protein_coding",
			HGVSp: "p.Gly20=", HGVSc: "c.60T>C",
		},
	}

	// Sort
	SortRows(rows)
	assert.Equal(t, int32(1), rows[0].ChromNumeric, "NRAS (chr1) should sort first")
	assert.Equal(t, int32(12), rows[1].ChromNumeric, "KRAS (chr12) should sort second")

	// Write to temp file
	f, err := os.CreateTemp(t.TempDir(), "test-*.parquet")
	require.NoError(t, err)

	w := NewWriter(f, 10) // small row group for test
	require.NoError(t, w.WriteRows(rows))
	require.NoError(t, w.Close())
	require.NoError(t, f.Close())

	// Read back
	rf, err := os.Open(f.Name())
	require.NoError(t, err)
	defer rf.Close()

	fi, err := rf.Stat()
	require.NoError(t, err)

	pf, err := pq.OpenFile(rf, fi.Size())
	require.NoError(t, err)

	assert.Equal(t, int64(2), pf.NumRows())

	// Read rows
	reader := pq.NewGenericReader[Row](pf)
	defer reader.Close()

	got := make([]Row, 2)
	n, err := reader.Read(got)
	require.NoError(t, err)
	assert.Equal(t, 2, n)

	// Verify NRAS row (should be first after sort)
	assert.Equal(t, "NRAS", got[0].GeneName)
	assert.Equal(t, int32(1), got[0].ChromNumeric)
	assert.Equal(t, int64(115256530), got[0].Pos)
	assert.Equal(t, "synonymous_variant", got[0].Consequence)

	// Verify KRAS row
	assert.Equal(t, "KRAS", got[1].GeneName)
	assert.Equal(t, int32(12), got[1].ChromNumeric)
	assert.Equal(t, int64(25245350), got[1].Pos)
	assert.Equal(t, "missense_variant", got[1].Consequence)
	assert.InDelta(t, float32(0.9876), got[1].AMScore, 0.001)
	assert.Equal(t, "ONCOGENE", got[1].OncokbGeneType)
	assert.True(t, got[1].IsCanonical)
}

func TestAnnotationToRow(t *testing.T) {
	ann := &annotate.Annotation{
		TranscriptID: "ENST00000311936.8",
		GeneName:     "KRAS",
		Consequence:  "missense_variant",
		Impact:       "MODERATE",
		IsCanonical:  true,
	}
	ann.SetExtra("oncokb", "gene_type", "ONCOGENE")

	row := AnnotationToRow("chrX", 12345, "A", "T", ann)
	assert.Equal(t, int32(23), row.ChromNumeric)
	assert.Equal(t, "chrX", row.Chrom)
	assert.Equal(t, int64(12345), row.Pos)
	assert.Equal(t, "KRAS", row.GeneName)
	assert.Equal(t, "ONCOGENE", row.OncokbGeneType)
	assert.True(t, row.IsCanonical)
}
