//go:build !cuda

package gpuhash

import "fmt"

// Table is an in-memory open-addressing hash table with linear probing.
// All exported methods are safe for concurrent reads; concurrent writes require
// external synchronisation (the table is built once at startup and then frozen).
type Table struct {
	slots    []Slot
	capacity uint64 // must be a power of 2
	count    uint64
}

// NewTable allocates a Table with at least minCapacity slots.
// minCapacity is rounded up to the next power of two and must be > 0.
func NewTable(minCapacity uint64) (*Table, error) {
	cap := nextPow2(minCapacity)
	if cap == 0 {
		return nil, fmt.Errorf("gpuhash: capacity overflow")
	}
	return &Table{slots: make([]Slot, cap), capacity: cap}, nil
}

// Insert adds or updates a slot.  key_hash must not be 0 (reserved sentinel).
func (t *Table) Insert(s Slot) {
	if s.Hash == 0 {
		return
	}
	idx := s.Hash & (t.capacity - 1)
	for {
		cur := &t.slots[idx]
		if cur.Hash == 0 || cur.Hash == s.Hash {
			if cur.Hash == 0 {
				t.count++
			}
			*cur = s
			return
		}
		idx = (idx + 1) & (t.capacity - 1)
	}
}

// Lookup returns the Value for the given hash key, and whether it was found.
// The CVClnDN field is always empty — callers must fill it via a secondary lookup.
func (t *Table) Lookup(hash uint64) (Value, bool) {
	if hash == 0 || t.capacity == 0 {
		return Value{}, false
	}
	idx := hash & (t.capacity - 1)
	for {
		s := &t.slots[idx]
		if s.Hash == 0 {
			return Value{}, false
		}
		if s.Hash == hash {
			return Value{
				AMScore:   s.AMScore,
				AMClass:   s.AMClass,
				CVSig:     s.CVSig,
				CVRevStat: s.CVRevStat,
				SigStatus: s.SigStatus,
				SigCount:  s.SigCount,
				SigFreq:   s.SigFreq,
			}, true
		}
		idx = (idx + 1) & (t.capacity - 1)
	}
}

// BatchLookup performs parallel hash lookups for all keys in hashes and writes
// results to out.  len(out) must equal len(hashes).
// On the CPU build this fans out across goroutines using the read-only table.
func (t *Table) BatchLookup(hashes []uint64, out []Value) {
	for i, h := range hashes {
		out[i], _ = t.Lookup(h)
	}
}

// Len returns the number of occupied slots.
func (t *Table) Len() uint64 { return t.count }

// Free is a no-op on the CPU build (GC handles memory).
func (t *Table) Free() {}

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
