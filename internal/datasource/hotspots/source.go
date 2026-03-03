package hotspots

import (
	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// Pre-built keys for Extra map (avoids string concatenation per annotation).
const (
	extraKeyHotspot = "hotspots.hotspot"
	extraKeyType    = "hotspots.type"
	extraKeyQValue  = "hotspots.qvalue"
)

// Source implements annotate.AnnotationSource for cancer hotspots.
type Source struct {
	store *Store
}

// NewSource creates an AnnotationSource backed by the given Store.
func NewSource(store *Store) *Source {
	return &Source{store: store}
}

func (s *Source) Name() string    { return "hotspots" }
func (s *Source) Version() string { return "v2" }

func (s *Source) Columns() []annotate.ColumnDef {
	return []annotate.ColumnDef{
		{Name: "hotspot", Description: "Y if position is a known cancer hotspot"},
		{Name: "type", Description: "Hotspot type: single residue, in-frame indel, 3d, splice"},
		{Name: "qvalue", Description: "Statistical significance (q-value)"},
	}
}

// Annotate marks annotations whose gene+protein position is a known hotspot.
func (s *Source) Annotate(_ *vcf.Variant, anns []*annotate.Annotation) {
	for _, ann := range anns {
		if ann.ProteinPosition == 0 || ann.GeneName == "" {
			continue
		}
		if h, ok := s.store.Lookup(ann.GeneName, ann.ProteinPosition); ok {
			ann.SetExtraKey(extraKeyHotspot, "Y")
			if h.Type != "" {
				ann.SetExtraKey(extraKeyType, h.Type)
			}
			ann.SetExtraKey(extraKeyQValue, FormatQValue(h.QValue))
		}
	}
}
