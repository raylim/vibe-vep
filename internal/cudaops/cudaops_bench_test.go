// Benchmarks comparing CPU vs GPU codon translation and reverse complement.
// Run without the cuda tag to benchmark the CPU fallback only:
//
//	go test ./internal/cudaops/ -bench=. -benchmem
//
// Run with the cuda tag to benchmark the GPU path (requires compiled libcudaops.so):
//
//	go test -tags cuda ./internal/cudaops/ -bench=. -benchmem

package cudaops_test

import (
	"testing"

	"github.com/inodb/vibe-vep/internal/cudaops"
)

// generateCodons builds a packed codon byte slice of n codons using the given bases.
func generateCodons(n int) []byte {
	bases := []byte("ACGT")
	out := make([]byte, 3*n)
	for i := range n {
		out[i*3+0] = bases[(i*7+0)%4]
		out[i*3+1] = bases[(i*7+1)%4]
		out[i*3+2] = bases[(i*7+2)%4]
	}
	return out
}

// generateSequence builds a random-ish DNA sequence of length n.
func generateSequence(n int) []byte {
	bases := []byte("ACGT")
	out := make([]byte, n)
	for i := range n {
		out[i] = bases[(i*13+7)%4]
	}
	return out
}

// ---- TranslateCodons benchmarks ----

func benchmarkTranslateCodons(b *testing.B, n int) {
	b.Helper()
	codons := generateCodons(n)
	b.ResetTimer()
	b.SetBytes(int64(3 * n))
	for b.Loop() {
		if _, err := cudaops.TranslateCodons(codons); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkTranslateCodons_1(b *testing.B)    { benchmarkTranslateCodons(b, 1) }
func BenchmarkTranslateCodons_64(b *testing.B)   { benchmarkTranslateCodons(b, 64) }
func BenchmarkTranslateCodons_1k(b *testing.B)   { benchmarkTranslateCodons(b, 1_000) }
func BenchmarkTranslateCodons_64k(b *testing.B)  { benchmarkTranslateCodons(b, 64_000) }
func BenchmarkTranslateCodons_1M(b *testing.B)   { benchmarkTranslateCodons(b, 1_000_000) }

// ---- TranslateCodonsPinned benchmarks (CUDA pinned / heap allocation) ----

func benchmarkTranslateCodonsPinned(b *testing.B, n int) {
	b.Helper()
	inBuf, err := cudaops.NewPinnedBuffer(3 * n)
	if err != nil {
		b.Fatal(err)
	}
	defer inBuf.Free()
	outBuf, err := cudaops.NewPinnedBuffer(n)
	if err != nil {
		b.Fatal(err)
	}
	defer outBuf.Free()

	// Fill input buffer with synthetic codons.
	bases := []byte("ACGT")
	in := inBuf.Bytes()
	for i := range n {
		in[i*3+0] = bases[(i*7+0)%4]
		in[i*3+1] = bases[(i*7+1)%4]
		in[i*3+2] = bases[(i*7+2)%4]
	}

	b.ResetTimer()
	b.SetBytes(int64(3 * n))
	for b.Loop() {
		if err := cudaops.TranslateCodonsPinned(inBuf, outBuf, n); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkTranslateCodonsPinned_1(b *testing.B)   { benchmarkTranslateCodonsPinned(b, 1) }
func BenchmarkTranslateCodonsPinned_64(b *testing.B)  { benchmarkTranslateCodonsPinned(b, 64) }
func BenchmarkTranslateCodonsPinned_1k(b *testing.B)  { benchmarkTranslateCodonsPinned(b, 1_000) }
func BenchmarkTranslateCodonsPinned_64k(b *testing.B) { benchmarkTranslateCodonsPinned(b, 64_000) }
func BenchmarkTranslateCodonsPinned_1M(b *testing.B)  { benchmarkTranslateCodonsPinned(b, 1_000_000) }

// ---- ReverseComplement benchmarks ----

func benchmarkReverseComplement(b *testing.B, n int) {
	b.Helper()
	seq := generateSequence(n)
	b.ResetTimer()
	b.SetBytes(int64(n))
	for b.Loop() {
		if _, err := cudaops.ReverseComplement(seq); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkReverseComplement_100bp(b *testing.B)    { benchmarkReverseComplement(b, 100) }
func BenchmarkReverseComplement_1kbp(b *testing.B)     { benchmarkReverseComplement(b, 1_000) }
func BenchmarkReverseComplement_100kbp(b *testing.B)   { benchmarkReverseComplement(b, 100_000) }
func BenchmarkReverseComplement_10Mbp(b *testing.B)    { benchmarkReverseComplement(b, 10_000_000) }

// ---- Correctness smoke tests ----

func TestTranslateCodons(t *testing.T) {
	codons := []byte("ATGGGG" + "TAAAAA") // ATG=M, GGG=G, TAA=*, AAA=K
	got, err := cudaops.TranslateCodons(codons)
	if err != nil {
		t.Fatal(err)
	}
	want := []byte{'M', 'G', '*', 'K'}
	if string(got) != string(want) {
		t.Errorf("TranslateCodons(%q) = %q, want %q", codons, got, want)
	}
}

func TestReverseComplement(t *testing.T) {
	cases := []struct{ in, want string }{
		{"ATGC", "GCAT"},
		{"AAAA", "TTTT"},
		{"GCGC", "GCGC"},
		{"A", "T"},
	}
	for _, tc := range cases {
		got, err := cudaops.ReverseComplement([]byte(tc.in))
		if err != nil {
			t.Fatalf("ReverseComplement(%q): %v", tc.in, err)
		}
		if string(got) != tc.want {
			t.Errorf("ReverseComplement(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}
