// Package parquet provides Parquet file export for variant annotations.
// The output files are sorted by (chrom_numeric, pos) to enable efficient
// row group pruning when queried via DuckDB-WASM in the browser.
package parquet

// Row represents a single variant annotation row in the Parquet file.
// Sort key: chrom_numeric (int32) + pos enables DuckDB row group stats pruning.
type Row struct {
	// Sort key
	ChromNumeric int32  `parquet:"chrom_numeric,zstd"`
	Pos          int64  `parquet:"pos,zstd"`
	Chrom        string `parquet:"chrom,dict,zstd"`
	Ref          string `parquet:"ref,zstd"`
	Alt          string `parquet:"alt,zstd"`

	// Transcript and gene
	TranscriptID string `parquet:"transcript_id,zstd"`
	GeneName     string `parquet:"gene_name,dict,zstd"`
	GeneID       string `parquet:"gene_id,zstd"`

	// Consequence
	Consequence     string `parquet:"consequence,dict,zstd"`
	Impact          string `parquet:"impact,dict,zstd"`
	CDSPosition     int64  `parquet:"cds_position,zstd"`
	ProteinPosition int64  `parquet:"protein_position,zstd"`
	AminoAcidChange string `parquet:"amino_acid_change,zstd"`
	CodonChange     string `parquet:"codon_change,zstd"`

	// Transcript details
	IsCanonicalMSK     bool `parquet:"is_canonical_msk,zstd"`
	IsCanonicalEnsembl bool `parquet:"is_canonical_ensembl,zstd"`
	Allele       string `parquet:"allele,zstd"`
	Biotype      string `parquet:"biotype,dict,zstd"`
	ExonNumber   string `parquet:"exon_number,zstd"`
	IntronNumber string `parquet:"intron_number,zstd"`
	CDNAPosition int64  `parquet:"cdna_position,zstd"`

	// HGVS notation
	HGVSp string `parquet:"hgvsp,zstd"`
	HGVSc string `parquet:"hgvsc,zstd"`

	// OncoKB
	OncokbGeneType string `parquet:"oncokb_gene_type,dict,zstd"`

	// AlphaMissense
	AMScore float32 `parquet:"am_score,zstd"`
	AMClass string  `parquet:"am_class,dict,zstd"`

	// ClinVar
	ClinvarClnSig    string `parquet:"clinvar_clnsig,dict,zstd"`
	ClinvarClnRevStat string `parquet:"clinvar_clnrevstat,dict,zstd"`
	ClinvarClnDN     string `parquet:"clinvar_clndn,zstd"`

	// Cancer Hotspots
	HotspotsHotspot string `parquet:"hotspots_hotspot,dict,zstd"`
	HotspotsType    string `parquet:"hotspots_type,dict,zstd"`
	HotspotsQValue  string `parquet:"hotspots_qvalue,zstd"`

	// SIGNAL
	SignalMutationStatus string `parquet:"signal_mutation_status,dict,zstd"`
	SignalCountCarriers  string `parquet:"signal_count_carriers,zstd"`
	SignalFrequency      string `parquet:"signal_frequency,zstd"`
}
