package signal

import (
	"strconv"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Pre-built keys for Extra map.
const (
	extraKeyMutationStatus = "signal.mutation_status"
	extraKeyCountCarriers  = "signal.count_carriers"
	extraKeyFrequency      = "signal.frequency"
)

// Source implements annotate.AnnotationSource for SIGNAL.
type Source struct {
	store *Store
}

// NewSource creates an AnnotationSource backed by the given Store.
func NewSource(store *Store) *Source {
	return &Source{store: store}
}

func (s *Source) Name() string    { return "signal" }
func (s *Source) Version() string { return "1.0" }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "mutation_status", Description: "Germline mutation status"},
		{Name: "count_carriers", Description: "Number of carriers in SIGNAL cohort"},
		{Name: "frequency", Description: "Overall allele frequency in SIGNAL cohort"},
	}
}

// Annotate adds SIGNAL germline mutation data to annotations.
func (s *Source) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {
	chrom := v.NormalizeChrom()
	entry, ok := s.store.Lookup(chrom, v.Pos, v.Ref, v.Alt)
	if !ok {
		return
	}
	for _, ann := range anns {
		ann.SetExtraKey(extraKeyMutationStatus, "germline")
		if entry.CountAll > 0 {
			ann.SetExtraKey(extraKeyCountCarriers, strconv.Itoa(entry.CountAll))
		}
		if entry.FreqAll > 0 {
			ann.SetExtraKey(extraKeyFrequency, FormatFreq(entry.FreqAll))
		}
	}
}
