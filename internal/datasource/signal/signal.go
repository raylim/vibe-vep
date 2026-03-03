// Package signal provides SIGNAL germline mutation frequency lookups.
// Data comes from the MSK SIGNAL database of germline mutations in cancer patients.
package signal

import (
	"bufio"
	"encoding/gob"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Entry represents a single SIGNAL variant annotation.
type Entry struct {
	Pos       int64
	Ref       string
	Alt       string
	Gene      string
	CountAll  int     // n_impact: total carrier count
	FreqAll   float64 // f_impact: overall frequency
	Biallelic float64 // f_biallelic: biallelic frequency
}

// Store holds SIGNAL data as sorted slices per chromosome.
type Store struct {
	data map[string][]Entry // chromosome (no "chr" prefix) → sorted entries
}

// Load parses a SIGNAL frequencies TSV file into sorted in-memory slices.
// Expected format: signaldb_all_variants_frequencies.txt
func Load(path string) (*Store, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open signal file: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)

	// Read header
	if !scanner.Scan() {
		return nil, fmt.Errorf("empty signal file")
	}
	header := strings.Split(scanner.Text(), "\t")
	col := indexColumns(header,
		"Hugo_Symbol", "Chromosome", "Start_Position",
		"Reference_Allele", "Alternate_Allele",
		"n_impact", "f_impact", "f_biallelic",
	)

	if col["Chromosome"] < 0 || col["Start_Position"] < 0 {
		return nil, fmt.Errorf("missing required columns: Chromosome, Start_Position")
	}

	data := make(map[string][]Entry)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		chromIdx := col["Chromosome"]
		posIdx := col["Start_Position"]
		if chromIdx >= len(fields) || posIdx >= len(fields) {
			continue
		}

		chrom := normalizeChrom(fields[chromIdx])
		pos, err := strconv.ParseInt(fields[posIdx], 10, 64)
		if err != nil {
			continue
		}

		entry := Entry{
			Pos:  pos,
			Ref:  getField(fields, col["Reference_Allele"]),
			Alt:  getField(fields, col["Alternate_Allele"]),
			Gene: getField(fields, col["Hugo_Symbol"]),
		}

		if idx := col["n_impact"]; idx >= 0 && idx < len(fields) {
			entry.CountAll, _ = strconv.Atoi(fields[idx])
		}
		if idx := col["f_impact"]; idx >= 0 && idx < len(fields) {
			entry.FreqAll, _ = strconv.ParseFloat(fields[idx], 64)
		}
		if idx := col["f_biallelic"]; idx >= 0 && idx < len(fields) {
			entry.Biallelic, _ = strconv.ParseFloat(fields[idx], 64)
		}

		// Convert SIGNAL convention: "-" means empty allele (insertion/deletion)
		if entry.Ref == "-" {
			entry.Ref = ""
		}
		if entry.Alt == "-" {
			entry.Alt = ""
		}

		data[chrom] = append(data[chrom], entry)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("read signal file: %w", err)
	}

	// Sort by position, then ref, then alt
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

// Lookup finds SIGNAL annotations for a specific variant.
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

// Count returns the total number of SIGNAL entries.
func (s *Store) Count() int {
	var n int
	for _, entries := range s.data {
		n += len(entries)
	}
	return n
}

func normalizeChrom(chrom string) string {
	return strings.TrimPrefix(chrom, "chr")
}

func getField(fields []string, idx int) string {
	if idx < 0 || idx >= len(fields) {
		return ""
	}
	return fields[idx]
}

// indexColumns maps column names to their indices in the header.
func indexColumns(header []string, names ...string) map[string]int {
	idx := make(map[string]int, len(names))
	for _, name := range names {
		idx[name] = -1
	}
	for i, col := range header {
		for _, name := range names {
			if col == name {
				idx[name] = i
			}
		}
	}
	return idx
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

// GobValid checks if a gob cache file exists and is newer than the source file.
func GobValid(gobPath, srcPath string) bool {
	gobInfo, err := os.Stat(gobPath)
	if err != nil {
		return false
	}
	srcInfo, err := os.Stat(srcPath)
	if err != nil {
		return false
	}
	return gobInfo.ModTime().After(srcInfo.ModTime())
}

// FormatFreq formats a frequency value for output.
func FormatFreq(v float64) string {
	if v == 0 {
		return "0"
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}
