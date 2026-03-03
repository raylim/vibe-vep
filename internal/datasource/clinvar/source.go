package clinvar

import (
	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Pre-built keys for Extra map.
const (
	extraKeyClnSig  = "clinvar.clnsig"
	extraKeyRevStat = "clinvar.clnrevstat"
	extraKeyClnDN   = "clinvar.clndn"
)

// Source implements annotate.AnnotationSource for ClinVar.
type Source struct {
	store *Store
}

// NewSource creates an AnnotationSource backed by the given Store.
func NewSource(store *Store) *Source {
	return &Source{store: store}
}

func (s *Source) Name() string    { return "clinvar" }
func (s *Source) Version() string { return "2024" }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "clnsig", Description: "Clinical significance (e.g. Pathogenic, Benign)"},
		{Name: "clnrevstat", Description: "Review status (e.g. reviewed_by_expert_panel)"},
		{Name: "clndn", Description: "Disease name(s)"},
	}
}

// Annotate adds ClinVar clinical significance to annotations.
func (s *Source) Annotate(v *vcf.Variant, anns []*annotate.Annotation) {
	chrom := v.NormalizeChrom()
	entry, ok := s.store.Lookup(chrom, v.Pos, v.Ref, v.Alt)
	if !ok {
		return
	}
	for _, ann := range anns {
		ann.SetExtraKey(extraKeyClnSig, entry.ClnSig)
		if entry.RevStat != "" {
			ann.SetExtraKey(extraKeyRevStat, entry.RevStat)
		}
		if entry.ClnDN != "" {
			ann.SetExtraKey(extraKeyClnDN, entry.ClnDN)
		}
	}
}
