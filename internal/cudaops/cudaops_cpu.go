// Package cudaops provides NVIDIA CUDA-accelerated genomic sequence operations.
// This file is a pure-Go stub used when the cuda build tag is NOT set.
// All operations fall back to CPU implementations.

//go:build !cuda

package cudaops

import "fmt"

// Init is a no-op in the CPU fallback build.
func Init() {}

// TranslateCodons translates codons using the CPU codon table.
// This fallback exists so code using cudaops compiles without CUDA.
func TranslateCodons(codons []byte) ([]byte, error) {
	if len(codons)%3 != 0 {
		return nil, fmt.Errorf("cudaops: codon slice length %d not divisible by 3", len(codons))
	}
	n := len(codons) / 3
	out := make([]byte, n)
	for i := range n {
		out[i] = translateCodonCPU(codons[i*3], codons[i*3+1], codons[i*3+2])
	}
	return out, nil
}

// PinnedBuffer is a regular heap buffer in the CPU fallback build.
type PinnedBuffer struct {
	data []byte
}

// NewPinnedBuffer allocates a regular heap buffer.
func NewPinnedBuffer(n int) (*PinnedBuffer, error) {
	return &PinnedBuffer{data: make([]byte, n)}, nil
}

// Free is a no-op for heap buffers.
func (b *PinnedBuffer) Free() {}

// Bytes returns the buffer as a byte slice.
func (b *PinnedBuffer) Bytes() []byte { return b.data }

// TranslateCodonsPinned translates using the CPU in the fallback build.
func TranslateCodonsPinned(in, out *PinnedBuffer, n int) error {
	if n <= 0 {
		return nil
	}
	if len(in.data) < 3*n || len(out.data) < n {
		return fmt.Errorf("cudaops: buffer too small")
	}
	for i := range n {
		out.data[i] = translateCodonCPU(in.data[i*3], in.data[i*3+1], in.data[i*3+2])
	}
	return nil
}

// ReverseComplement reverse-complements a DNA sequence on the CPU.
func ReverseComplement(seq []byte) ([]byte, error) {
	n := len(seq)
	out := make([]byte, n)
	for i := range n {
		out[i] = complementCPU(seq[n-1-i])
	}
	return out, nil
}

var codonTableCPU = [64]byte{
	// Indices 0–15: first base A (A=0)
	'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
	// Indices 16–31: first base C
	'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
	// Indices 32–47: first base G
	'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
	// Indices 48–63: first base T
	'*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F',
}

func base2idx(b byte) byte {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	default:
		return 0
	}
}

func translateCodonCPU(b0, b1, b2 byte) byte {
	idx := (base2idx(b0) << 4) | (base2idx(b1) << 2) | base2idx(b2)
	return codonTableCPU[idx]
}

func complementCPU(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'G':
		return 'C'
	case 'C':
		return 'G'
	case 'a':
		return 't'
	case 't':
		return 'a'
	case 'g':
		return 'c'
	case 'c':
		return 'g'
	default:
		return 'N'
	}
}
