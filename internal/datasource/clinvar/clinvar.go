// Package clinvar provides ClinVar clinical significance lookups.
// Data comes from the ClinVar VCF file (NCBI).
package clinvar

import (
	"bufio"
	"compress/gzip"
	"encoding/gob"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Entry represents a single ClinVar variant annotation.
type Entry struct {
	Pos    int64
	Ref    string
	Alt    string
	ClnSig string // Clinical significance (e.g., "Pathogenic")
	RevStat string // Review status
	ClnDN  string // Disease name(s)
}

// Store holds ClinVar data as sorted slices per chromosome.
type Store struct {
	data map[string][]Entry // chromosome (no "chr" prefix) → sorted entries
}

// Load parses a ClinVar VCF file (gzipped or plain) into sorted in-memory slices.
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open clinvar file: %w", err)
	}
	defer f.Close()

	var scanner *bufio.Scanner
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("open gzip reader: %w", err)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	} else {
		scanner = bufio.NewScanner(f)
	}
	scanner.Buffer(make([]byte, 4*1024*1024), 4*1024*1024)

	data := make(map[string][]Entry)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}

		entry, chrom, ok := parseVCFLine(line)
		if !ok {
			continue
		}
		data[chrom] = append(data[chrom], entry)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("read clinvar file: %w", err)
	}

	// Sort each chromosome's entries by position
	for chrom, entries := range data {
		sort.Slice(entries, func(i, j int) bool {
			if entries[i].Pos != entries[j].Pos {
				return entries[i].Pos < entries[j].Pos
			}
			if entries[i].Ref != entries[j].Ref {
				return entries[i].Ref < entries[j].Ref
			}
			return entries[i].Alt < entries[j].Alt
		})
		data[chrom] = entries
	}

	return &Store{data: data}, nil
}

// parseVCFLine parses a single VCF data line into a ClinVar Entry.
// Returns the entry, normalized chromosome, and whether parsing succeeded.
func parseVCFLine(line string) (Entry, string, bool) {
	// VCF: CHROM POS ID REF ALT QUAL FILTER INFO ...
	fields := strings.SplitN(line, "\t", 9)
	if len(fields) < 8 {
		return Entry{}, "", false
	}

	chrom := normalizeChrom(fields[0])
	pos, err := strconv.ParseInt(fields[1], 10, 64)
	if err != nil {
		return Entry{}, "", false
	}

	ref := fields[3]
	alt := fields[4]
	info := fields[7]

	// Handle multi-allelic: take first ALT only
	if idx := strings.IndexByte(alt, ','); idx >= 0 {
		alt = alt[:idx]
	}

	entry := Entry{
		Pos:     pos,
		Ref:     ref,
		Alt:     alt,
		ClnSig:  extractInfo(info, "CLNSIG="),
		RevStat: extractInfo(info, "CLNREVSTAT="),
		ClnDN:   extractInfo(info, "CLNDN="),
	}

	// Skip entries without clinical significance
	if entry.ClnSig == "" {
		return Entry{}, "", false
	}

	return entry, chrom, true
}

// extractInfo extracts a value from a VCF INFO field.
func extractInfo(info, key string) string {
	idx := strings.Index(info, key)
	if idx < 0 {
		return ""
	}
	val := info[idx+len(key):]
	if end := strings.IndexByte(val, ';'); end >= 0 {
		val = val[:end]
	}
	return val
}

// Lookup finds ClinVar annotations for a specific variant.
func (s *Store) Lookup(chrom string, pos int64, ref, alt string) (Entry, bool) {
	chrom = normalizeChrom(chrom)
	entries := s.data[chrom]
	if len(entries) == 0 {
		return Entry{}, false
	}

	// Binary search for position
	i := sort.Search(len(entries), func(i int) bool {
		return entries[i].Pos >= pos
	})

	// Scan forward through entries at this position
	for ; i < len(entries) && entries[i].Pos == pos; i++ {
		if entries[i].Ref == ref && entries[i].Alt == alt {
			return entries[i], true
		}
	}
	return Entry{}, false
}

// Count returns the total number of ClinVar entries.
func (s *Store) Count() int {
	var n int
	for _, entries := range s.data {
		n += len(entries)
	}
	return n
}

// normalizeChrom removes "chr" prefix for consistent lookups.
func normalizeChrom(chrom string) string {
	return strings.TrimPrefix(chrom, "chr")
}

// gobData is the serialization format for gob cache.
type gobData struct {
	Data map[string][]Entry
}

// SaveGob serializes the store to a gob file for fast reload.
func (s *Store) SaveGob(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create gob file: %w", err)
	}

	if err := gob.NewEncoder(f).Encode(gobData{Data: s.data}); err != nil {
		f.Close()
		os.Remove(path)
		return fmt.Errorf("encode gob: %w", err)
	}
	return f.Close()
}

// LoadGob deserializes a store from a gob file.
func LoadGob(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open gob file: %w", err)
	}
	defer f.Close()

	var d gobData
	if err := gob.NewDecoder(f).Decode(&d); err != nil {
		return nil, fmt.Errorf("decode gob: %w", err)
	}
	return &Store{data: d.Data}, nil
}

// GobValid checks if a gob cache file exists and is newer than the source VCF.
func GobValid(gobPath, vcfPath string) bool {
	gobInfo, err := os.Stat(gobPath)
	if err != nil {
		return false
	}
	vcfInfo, err := os.Stat(vcfPath)
	if err != nil {
		return false
	}
	return gobInfo.ModTime().After(vcfInfo.ModTime())
}
