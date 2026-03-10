// Package cudaops provides NVIDIA CUDA-accelerated genomic sequence operations.
//
// Build with the cuda tag to enable GPU acceleration:
//
//	CGO_ENABLED=1 CGO_CFLAGS="-I/usr/local/cuda/include" \
//	  CGO_LDFLAGS="-L/usr/local/cuda/lib64 -lcudart" \
//	  go build -tags cuda ./...
//
// The shared library must be compiled first:
//
//	nvcc -O3 -arch=sm_80 --compiler-options -fPIC \
//	  -shared -o internal/cudaops/libcudaops.so \
//	  internal/cudaops/codon_cuda.cu
//
// At runtime the library is loaded from the directory containing the binary
// or from LD_LIBRARY_PATH.

//go:build cuda

package cudaops

/*
#cgo CFLAGS:  -I/usr/local/cuda/include
#cgo LDFLAGS: -L${SRCDIR} -lcudaops -Wl,-rpath,${SRCDIR} -L/usr/local/cuda/lib64 -lcudart

#include <stdlib.h>

void  cudaops_init_complement_table();
int   cudaops_translate_codons(const char* in, char* out, int n);
int   cudaops_translate_codons_pinned(const char* in, char* out, int n);
int   cudaops_reverse_complement(const char* in, char* out, int n);
void* cudaops_alloc_pinned(int n);
void  cudaops_free_pinned(void* ptr);
*/
import "C"
import (
	"fmt"
	"sync"
	"unsafe"
)

var initOnce sync.Once

// Init initialises the CUDA constant memory tables.
// It is called automatically by the first operation, but may be called
// explicitly to front-load the initialisation cost.
func Init() {
	initOnce.Do(func() {
		C.cudaops_init_complement_table()
	})
}

// TranslateCodons translates a flat byte slice of packed codons (3 bytes each)
// to amino acid single-letter codes using the GPU.
// len(codons) must be divisible by 3; the returned slice has len/3 elements.
func TranslateCodons(codons []byte) ([]byte, error) {
	Init()
	if len(codons) == 0 {
		return nil, nil
	}
	if len(codons)%3 != 0 {
		return nil, fmt.Errorf("cudaops: codon slice length %d not divisible by 3", len(codons))
	}
	n := len(codons) / 3
	out := make([]byte, n)

	rc := C.cudaops_translate_codons(
		(*C.char)(unsafe.Pointer(&codons[0])),
		(*C.char)(unsafe.Pointer(&out[0])),
		C.int(n),
	)
	if rc != 0 {
		return nil, fmt.Errorf("cudaops: translate_codons CUDA error %d", rc)
	}
	return out, nil
}

// PinnedBuffer is a CUDA pinned-memory buffer for zero-copy transfers.
// Use it to amortise PCIe overhead when calling TranslateCodonsPinned repeatedly.
type PinnedBuffer struct {
	ptr  unsafe.Pointer
	size int
}

// NewPinnedBuffer allocates n bytes of pinned host memory.
func NewPinnedBuffer(n int) (*PinnedBuffer, error) {
	p := C.cudaops_alloc_pinned(C.int(n))
	if p == nil {
		return nil, fmt.Errorf("cudaops: failed to allocate %d bytes of pinned memory", n)
	}
	return &PinnedBuffer{ptr: p, size: n}, nil
}

// Free releases the pinned memory.
func (b *PinnedBuffer) Free() {
	if b.ptr != nil {
		C.cudaops_free_pinned(b.ptr)
		b.ptr = nil
	}
}

// Bytes returns the buffer as a byte slice.
func (b *PinnedBuffer) Bytes() []byte {
	return unsafe.Slice((*byte)(b.ptr), b.size)
}

// TranslateCodonsPinned is like TranslateCodons but uses caller-provided
// pinned input/output buffers to minimise PCIe transfer latency.
// in and out must be PinnedBuffers of size 3*n and n respectively.
func TranslateCodonsPinned(in, out *PinnedBuffer, n int) error {
	Init()
	if n <= 0 {
		return nil
	}
	if in.size < 3*n || out.size < n {
		return fmt.Errorf("cudaops: buffer too small (in=%d, out=%d, n=%d)", in.size, out.size, n)
	}
	rc := C.cudaops_translate_codons_pinned(
		(*C.char)(in.ptr),
		(*C.char)(out.ptr),
		C.int(n),
	)
	if rc != 0 {
		return fmt.Errorf("cudaops: translate_codons_pinned CUDA error %d", rc)
	}
	return nil
}

// ReverseComplement reverse-complements a DNA sequence using the GPU.
// Returns a new byte slice of the same length.
func ReverseComplement(seq []byte) ([]byte, error) {
	Init()
	if len(seq) == 0 {
		return nil, nil
	}
	out := make([]byte, len(seq))
	rc := C.cudaops_reverse_complement(
		(*C.char)(unsafe.Pointer(&seq[0])),
		(*C.char)(unsafe.Pointer(&out[0])),
		C.int(len(seq)),
	)
	if rc != 0 {
		return nil, fmt.Errorf("cudaops: reverse_complement CUDA error %d", rc)
	}
	return out, nil
}
