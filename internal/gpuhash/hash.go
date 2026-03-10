package gpuhash

import (
	"encoding/binary"
)

const (
	fnvOffset uint64 = 14695981039346656037
	fnvPrime  uint64 = 1099511628211
)

// HashKey computes the FNV-1a 64-bit hash of the canonical annotation lookup
// key: chrom\x00 + pos(4-byte little-endian) + \x00 + ref + \x00 + alt.
//
// Chromosome is expected without the "chr" prefix (e.g. "1", "X").
// Hash values 0 is reserved as the empty-slot sentinel; the function maps 0 to
// the constant emptyBias so callers never need to special-case it.
func HashKey(chrom string, pos int64, ref, alt string) uint64 {
	h := fnvOffset
	for i := 0; i < len(chrom); i++ {
		h ^= uint64(chrom[i])
		h *= fnvPrime
	}
	h ^= 0
	h *= fnvPrime
	var pb [4]byte
	binary.LittleEndian.PutUint32(pb[:], uint32(pos))
	for _, b := range pb {
		h ^= uint64(b)
		h *= fnvPrime
	}
	h ^= 0
	h *= fnvPrime
	for i := 0; i < len(ref); i++ {
		h ^= uint64(ref[i])
		h *= fnvPrime
	}
	h ^= 0
	h *= fnvPrime
	for i := 0; i < len(alt); i++ {
		h ^= uint64(alt[i])
		h *= fnvPrime
	}
	if h == 0 {
		return fnvPrime // bias away from reserved sentinel
	}
	return h
}
