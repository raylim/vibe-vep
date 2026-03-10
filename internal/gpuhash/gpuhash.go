//go:build cuda

// Package gpuhash provides GPU-accelerated batch hash table lookups.
// Build with -tags cuda to enable the CUDA path; otherwise the pure-Go
// in-memory table in gpuhash_cpu.go is used transparently.
//
// The shared library must be compiled before building with -tags cuda:
//
//	make gpuhash
package gpuhash

/*
#cgo CFLAGS: -I/usr/local/cuda-12.3/include
#cgo LDFLAGS: -L${SRCDIR} -lgpuhash -Wl,-rpath,${SRCDIR} -L/usr/local/cuda-12.3/lib64 -lcudart

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct __attribute__((packed, aligned(8))) {
    uint64_t hash;
    float    am_score;
    uint8_t  am_class;
    uint8_t  cv_sig;
    uint8_t  cv_revstat;
    uint8_t  sig_status;
    uint32_t sig_count;
    float    sig_freq;
    uint8_t  reserved[8];
} CSlot;

typedef struct __attribute__((packed, aligned(4))) {
    float    am_score;
    uint8_t  am_class;
    uint8_t  cv_sig;
    uint8_t  cv_revstat;
    uint8_t  sig_status;
    uint32_t sig_count;
    float    sig_freq;
    uint8_t  found;
    uint8_t  pad[3];
} CResult;

extern void   *gpuhash_create(uint64_t capacity);
extern int     gpuhash_load(void *d_table, const CSlot *h_slots, uint64_t count, uint64_t offset);
extern void    gpuhash_destroy(void *d_table);
extern int     gpuhash_batch_lookup(void *d_table, uint64_t capacity, void *d_hashes, void *d_results, uint32_t n);
extern int     gpuhash_sync(void);
extern void   *gpuhash_alloc_device(size_t bytes);
extern void    gpuhash_free_device(void *p);
extern int     gpuhash_copy_hashes_to_device(void *d_buf, const uint64_t *h_hashes, uint32_t n);
extern int     gpuhash_copy_results_to_host(void *h_buf, const void *d_buf, size_t result_bytes);
extern size_t  gpuhash_result_size(void);
*/
import "C"

import (
	"fmt"
	"unsafe"
)

// Table wraps a GPU-resident open-addressing hash table.
type Table struct {
	dTable   unsafe.Pointer // device memory for Slot array
	dHashes  unsafe.Pointer // device buffer for query hashes (reused across calls)
	dResults unsafe.Pointer // device buffer for results (reused across calls)
	capacity uint64
	count    uint64
	bufSize  uint32 // current size of the device hash/result buffers (in entries)
}

// NewTable allocates a GPU hash table with at least minCapacity slots.
// minCapacity is rounded up to the next power of two.
func NewTable(minCapacity uint64) (*Table, error) {
	cap := nextPow2(minCapacity)
	if cap == 0 {
		return nil, fmt.Errorf("gpuhash: capacity overflow")
	}
	d := C.gpuhash_create(C.uint64_t(cap))
	if d == nil {
		return nil, fmt.Errorf("gpuhash: cudaMalloc failed for %d slots", cap)
	}
	return &Table{dTable: d, capacity: cap}, nil
}

// LoadSlots uploads a slice of Slots into the device table starting at offset.
// Call this after NewTable to populate the table from host-side data.
func (t *Table) LoadSlots(slots []Slot, offset uint64) error {
	if len(slots) == 0 {
		return nil
	}
	rc := C.gpuhash_load(
		t.dTable,
		(*C.CSlot)(unsafe.Pointer(&slots[0])),
		C.uint64_t(len(slots)),
		C.uint64_t(offset),
	)
	if rc != 0 {
		return fmt.Errorf("gpuhash: cudaMemcpy to device failed")
	}
	t.count += uint64(len(slots))
	return nil
}

// BatchLookup performs parallel GPU hash lookups for all hashes and writes
// results to out.  len(out) must equal len(hashes).
// Device buffers are reused/grown as needed.
func (t *Table) BatchLookup(hashes []uint64, out []Value) error {
	n := uint32(len(hashes))
	if n == 0 {
		return nil
	}

	if err := t.ensureDeviceBuffers(n); err != nil {
		return err
	}

	// Copy hashes to device.
	if rc := C.gpuhash_copy_hashes_to_device(
		t.dHashes,
		(*C.uint64_t)(unsafe.Pointer(&hashes[0])),
		C.uint32_t(n),
	); rc != 0 {
		return fmt.Errorf("gpuhash: copy hashes to device failed")
	}

	// Launch kernel.
	if rc := C.gpuhash_batch_lookup(
		t.dTable,
		C.uint64_t(t.capacity),
		t.dHashes,
		t.dResults,
		C.uint32_t(n),
	); rc != 0 {
		return fmt.Errorf("gpuhash: kernel launch failed")
	}

	// Copy results back to host.
	resultBytes := uintptr(n) * uintptr(C.gpuhash_result_size())
	hostBuf := make([]C.CResult, n)
	if rc := C.gpuhash_copy_results_to_host(
		unsafe.Pointer(&hostBuf[0]),
		t.dResults,
		C.size_t(resultBytes),
	); rc != 0 {
		return fmt.Errorf("gpuhash: copy results to host failed")
	}

	if rc := C.gpuhash_sync(); rc != 0 {
		return fmt.Errorf("gpuhash: device synchronize failed")
	}

	// Decode C results to Go Values.
	for i := uint32(0); i < n; i++ {
		r := &hostBuf[i]
		if r.found == 0 {
			out[i] = Value{}
			continue
		}
		out[i] = Value{
			AMScore:   float32(r.am_score),
			AMClass:   AMClass(r.am_class),
			CVSig:     CVSig(r.cv_sig),
			CVRevStat: CVRevStat(r.cv_revstat),
			SigStatus: SigStatus(r.sig_status),
			SigCount:  uint32(r.sig_count),
			SigFreq:   float32(r.sig_freq),
		}
	}
	return nil
}

// Len returns the number of occupied slots.
func (t *Table) Len() uint64 { return t.count }

// Free releases all GPU memory held by this table.
func (t *Table) Free() {
	if t.dTable != nil {
		C.gpuhash_destroy(t.dTable)
		t.dTable = nil
	}
	if t.dHashes != nil {
		C.gpuhash_free_device(t.dHashes)
		t.dHashes = nil
	}
	if t.dResults != nil {
		C.gpuhash_free_device(t.dResults)
		t.dResults = nil
	}
}

// ensureDeviceBuffers grows the device hash/result buffers to fit n entries.
func (t *Table) ensureDeviceBuffers(n uint32) error {
	if n <= t.bufSize {
		return nil
	}
	if t.dHashes != nil {
		C.gpuhash_free_device(t.dHashes)
	}
	if t.dResults != nil {
		C.gpuhash_free_device(t.dResults)
	}
	hashBytes := C.size_t(uintptr(n) * 8)
	resultBytes := C.size_t(uintptr(n) * uintptr(C.gpuhash_result_size()))
	t.dHashes = C.gpuhash_alloc_device(hashBytes)
	t.dResults = C.gpuhash_alloc_device(resultBytes)
	if t.dHashes == nil || t.dResults == nil {
		return fmt.Errorf("gpuhash: device buffer allocation failed for %d entries", n)
	}
	t.bufSize = n
	return nil
}

// nextPow2 returns the smallest power of two >= n.
func nextPow2(n uint64) uint64 {
	if n == 0 {
		return 1
	}
	n--
	n |= n >> 1
	n |= n >> 2
	n |= n >> 4
	n |= n >> 8
	n |= n >> 16
	n |= n >> 32
	return n + 1
}
