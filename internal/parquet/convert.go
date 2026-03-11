package parquet

import (
	"fmt"
	"strings"

	"github.com/inodb/vibe-vep/internal/annotate"
)

// ChromToNumeric converts a chromosome name to a sortable integer.
// 1-22 map to 1-22, X=23, Y=24, MT/M=25. Unknown chromosomes return 0.
func ChromToNumeric(chrom string) int32 {
	// Strip "chr" prefix
	c := chrom
	if len(c) > 3 && strings.EqualFold(c[:3], "chr") {
		c = c[3:]
	}

	switch strings.ToUpper(c) {
	case "X":
		return 23
	case "Y":
		return 24
	case "MT", "M":
		return 25
	}

	var n int32
	if _, err := fmt.Sscanf(c, "%d", &n); err == nil && n >= 1 && n <= 22 {
		return n
	}
	return 0
}

// AnnotationToRow converts a variant's annotation to a Parquet Row.
// chrom should be the normalized chromosome (without "chr" prefix).
func AnnotationToRow(chrom string, pos int64, ref, alt string, ann *annotate.Annotation) Row {
	var amScore float32
	if s := ann.GetExtra("alphamissense", "score"); s != "" {
		fmt.Sscanf(s, "%f", &amScore)
	}

	return Row{
		ChromNumeric: ChromToNumeric(chrom),
		Pos:          pos,
		Chrom:        chrom,
		Ref:          ref,
		Alt:          alt,

		TranscriptID: ann.TranscriptID,
		GeneName:     ann.GeneName,
		GeneID:       ann.GeneID,

		Consequence:     ann.Consequence,
		Impact:          ann.Impact,
		CDSPosition:     ann.CDSPosition,
		ProteinPosition: ann.ProteinPosition,
		AminoAcidChange: ann.AminoAcidChange,
		CodonChange:     ann.CodonChange,

		IsCanonicalMSK:     ann.IsCanonicalMSK,
		IsCanonicalEnsembl: ann.IsCanonicalEnsembl,
		Allele:       ann.Allele,
		Biotype:      ann.Biotype,
		ExonNumber:   ann.ExonNumber,
		IntronNumber: ann.IntronNumber,
		CDNAPosition: ann.CDNAPosition,

		HGVSp: ann.HGVSp,
		HGVSc: ann.HGVSc,

		OncokbGeneType: ann.GetExtra("oncokb", "gene_type"),

		AMScore: amScore,
		AMClass: ann.GetExtra("alphamissense", "class"),

		ClinvarClnSig:    ann.GetExtra("clinvar", "clnsig"),
		ClinvarClnRevStat: ann.GetExtra("clinvar", "clnrevstat"),
		ClinvarClnDN:     ann.GetExtra("clinvar", "clndn"),

		HotspotsHotspot: ann.GetExtra("hotspots", "hotspot"),
		HotspotsType:    ann.GetExtra("hotspots", "type"),
		HotspotsQValue:  ann.GetExtra("hotspots", "qvalue"),

		SignalMutationStatus: ann.GetExtra("signal", "mutation_status"),
		SignalCountCarriers:  ann.GetExtra("signal", "count_carriers"),
		SignalFrequency:      ann.GetExtra("signal", "frequency"),
	}
}
