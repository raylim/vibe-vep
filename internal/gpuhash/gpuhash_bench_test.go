package gpuhash

import (
	"math/rand"
	"testing"
)

// BenchmarkHashKey measures the FNV-1a hash of a typical SNV key.
func BenchmarkHashKey(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = HashKey("17", 7674220, "G", "A")
	}
}

// BenchmarkTableLookupHit benchmarks a single successful hash table lookup.
func BenchmarkTableLookupHit(b *testing.B) {
	tbl := mustNewTable(b, 1024)
	defer tbl.Free()
	h := insertRandom(tbl, 512)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = tbl.Lookup(h)
	}
}

// BenchmarkBatchLookup_1k measures BatchLookup with 1 000 queries.
func BenchmarkBatchLookup_1k(b *testing.B) {
	benchBatch(b, 1_000)
}

// BenchmarkBatchLookup_64k measures BatchLookup with 64 000 queries.
func BenchmarkBatchLookup_64k(b *testing.B) {
	benchBatch(b, 64_000)
}

// BenchmarkBatchLookup_1M measures BatchLookup with 1 000 000 queries.
func BenchmarkBatchLookup_1M(b *testing.B) {
	benchBatch(b, 1_000_000)
}

func benchBatch(b *testing.B, n int) {
	b.Helper()
	tbl := mustNewTable(b, uint64(n*2))
	defer tbl.Free()

	hashes := make([]uint64, n)
	for i := range hashes {
		h := uint64(rand.Int63())
		if h == 0 {
			h = 1
		}
		tbl.Insert(Slot{Hash: h, AMScore: float32(i) / float32(n), AMClass: AMClassAmbiguous})
		hashes[i] = h
	}

	results := make([]Value, n)
	b.ResetTimer()
	b.ReportAllocs()
	b.SetBytes(int64(n * 8)) // bytes of hash input
	for i := 0; i < b.N; i++ {
		tbl.BatchLookup(hashes, results)
	}
}

func mustNewTable(tb testing.TB, cap uint64) *Table {
	tb.Helper()
	tbl, err := NewTable(cap)
	if err != nil {
		tb.Fatalf("NewTable: %v", err)
	}
	return tbl
}

// insertRandom inserts n random slots and returns the last hash inserted.
func insertRandom(tbl *Table, n int) uint64 {
	var last uint64
	for i := 0; i < n; i++ {
		h := uint64(rand.Int63())
		if h == 0 {
			h = 1
		}
		tbl.Insert(Slot{Hash: h, AMScore: float32(i)})
		last = h
	}
	return last
}
