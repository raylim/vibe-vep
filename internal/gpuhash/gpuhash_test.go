package gpuhash

import (
	"testing"
)

func TestHashKey(t *testing.T) {
	// Same inputs must produce same hash.
	h1 := HashKey("1", 12345, "A", "G")
	h2 := HashKey("1", 12345, "A", "G")
	if h1 != h2 {
		t.Fatalf("HashKey not deterministic: got %d and %d", h1, h2)
	}

	// Different inputs must (almost certainly) produce different hashes.
	h3 := HashKey("1", 12345, "A", "T")
	if h1 == h3 {
		t.Fatalf("HashKey collision between A>G and A>T at same position")
	}

	// Hash must never be 0 (reserved sentinel).
	tests := [][4]string{
		{"1", "0", "", ""},
		{"X", "999999999", "ACGT", "T"},
		{"MT", "1", "A", "G"},
	}
	for _, tc := range tests {
		pos := int64(0)
		for _, c := range tc[1] {
			pos = pos*10 + int64(c-'0')
		}
		h := HashKey(tc[0], pos, tc[2], tc[3])
		if h == 0 {
			t.Errorf("HashKey returned 0 for %v", tc)
		}
	}
}

func TestTableInsertLookup(t *testing.T) {
	tbl, err := NewTable(16)
	if err != nil {
		t.Fatal(err)
	}
	defer tbl.Free()

	// Insert a slot.
	key := HashKey("7", 140453136, "A", "T") // BRAF V600E-ish
	s := Slot{
		Hash:    key,
		AMScore: 0.93,
		AMClass: AMClassLikelyPathogenic,
		CVSig:   CVSigPathogenic,
	}
	tbl.Insert(s)

	v, ok := tbl.Lookup(key)
	if !ok {
		t.Fatal("lookup returned not found after insert")
	}
	if v.AMScore != 0.93 {
		t.Errorf("AMScore: got %v, want 0.93", v.AMScore)
	}
	if v.AMClass != AMClassLikelyPathogenic {
		t.Errorf("AMClass: got %v, want LikelyPathogenic", v.AMClass)
	}
	if v.CVSig != CVSigPathogenic {
		t.Errorf("CVSig: got %v, want Pathogenic", v.CVSig)
	}

	// Non-existent key returns not found.
	missing := HashKey("7", 999999, "C", "T")
	if _, ok := tbl.Lookup(missing); ok {
		t.Error("lookup returned found for a key that was not inserted")
	}
}

func TestTableBatchLookup(t *testing.T) {
	tbl, err := NewTable(64)
	if err != nil {
		t.Fatal(err)
	}
	defer tbl.Free()

	type entry struct {
		chrom string
		pos   int64
		ref   string
		alt   string
		score float32
	}
	entries := []entry{
		{"1", 100, "A", "G", 0.1},
		{"2", 200, "C", "T", 0.5},
		{"3", 300, "G", "A", 0.9},
	}
	for _, e := range entries {
		h := HashKey(e.chrom, e.pos, e.ref, e.alt)
		tbl.Insert(Slot{Hash: h, AMScore: e.score, AMClass: AMClassAmbiguous})
	}

	hashes := make([]uint64, 4)
	hashes[0] = HashKey("1", 100, "A", "G")
	hashes[1] = HashKey("99", 999, "T", "C") // not inserted
	hashes[2] = HashKey("2", 200, "C", "T")
	hashes[3] = HashKey("3", 300, "G", "A")

	results := make([]Value, 4)
	tbl.BatchLookup(hashes, results)

	if results[0].AMScore != 0.1 {
		t.Errorf("[0] AMScore: got %v, want 0.1", results[0].AMScore)
	}
	if results[1].HasAM() {
		t.Error("[1] expected miss, got hit")
	}
	if results[2].AMScore != 0.5 {
		t.Errorf("[2] AMScore: got %v, want 0.5", results[2].AMScore)
	}
	if results[3].AMScore != 0.9 {
		t.Errorf("[3] AMScore: got %v, want 0.9", results[3].AMScore)
	}
}

func TestEncodeDecodeAMClass(t *testing.T) {
	cases := []struct {
		s string
		c AMClass
	}{
		{"likely_benign", AMClassLikelyBenign},
		{"ambiguous", AMClassAmbiguous},
		{"likely_pathogenic", AMClassLikelyPathogenic},
		{"", AMClassAbsent},
		{"unknown_value", AMClassAbsent},
	}
	for _, tc := range cases {
		got := EncodeAMClass(tc.s)
		if got != tc.c {
			t.Errorf("EncodeAMClass(%q): got %d, want %d", tc.s, got, tc.c)
		}
		if tc.c != AMClassAbsent {
			if AMClassString[tc.c] != tc.s {
				t.Errorf("AMClassString[%d]: got %q, want %q", tc.c, AMClassString[tc.c], tc.s)
			}
		}
	}
}

func TestEncodeCVSig(t *testing.T) {
	cases := []string{
		"Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign",
		"Uncertain_significance", "Conflicting_interpretations_of_pathogenicity",
	}
	for _, s := range cases {
		code := EncodeCVSig(s)
		if code == CVSigAbsent {
			t.Errorf("EncodeCVSig(%q) returned Absent unexpectedly", s)
		}
		decoded, ok := CVSigStrings[code]
		if !ok {
			t.Errorf("CVSigStrings missing key for code %d (input %q)", code, s)
			continue
		}
		if decoded != s {
			t.Errorf("round-trip %q → %d → %q", s, code, decoded)
		}
	}
}
