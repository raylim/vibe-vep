// Package mane provides MANE Select transcript mapping between RefSeq and Ensembl.
//
// MANE (Matched Annotation from NCBI and EMBL-EBI) Select transcripts represent
// a consensus between RefSeq and Ensembl for every protein-coding gene. For MANE
// Select transcripts, both NM_ and ENST accessions encode identical protein
// sequences, so HGVSp comparisons across RefSeq-based and Ensembl-based tools are
// directly meaningful.
package mane

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// Map provides MANE Select transcript lookups by RefSeq accession.
//
// It maps both unversioned accessions (e.g. "NM_000546") to their Ensembl
// counterparts and retains the versioned RefSeq accession (e.g. "NM_000546.6")
// from the MANE summary file for exact-version comparisons.
type Map struct {
	toENST      map[string]string // unversioned NM_ → unversioned ENST
	nmVersioned map[string]string // unversioned NM_ → versioned NM_ from MANE file
}

// ENST returns the unversioned Ensembl transcript for a RefSeq accession.
// The refseq argument may be versioned (e.g. "NM_000546.6") or not.
// Returns ("", false) if not present.
func (m Map) ENST(refseq string) (string, bool) {
	k := stripVersion(refseq)
	v, ok := m.toENST[k]
	return v, ok
}

// HasRefSeq reports whether the RefSeq accession (versioned or not) is in MANE.
func (m Map) HasRefSeq(refseq string) bool {
	_, ok := m.toENST[stripVersion(refseq)]
	return ok
}

// HasExactVersion reports whether refseq is a MANE Select transcript and its
// version exactly matches the version recorded in the MANE summary file.
// For example, "NM_000546.6" returns true only if MANE lists NM_000546.6.
// An unversioned accession always returns false.
func (m Map) HasExactVersion(refseq string) bool {
	maneVersioned, ok := m.nmVersioned[stripVersion(refseq)]
	if !ok {
		return false
	}
	return refseq == maneVersioned
}

// Len returns the number of MANE Select transcripts in the map.
func (m Map) Len() int {
	return len(m.toENST)
}

// Load reads a MANE summary file (plain or gzipped) and returns the mapping.
// The file is the NCBI MANE summary (e.g. MANE.GRCh38.v1.5.summary.txt.gz).
func Load(path string) (Map, error) {
	f, err := os.Open(path)
	if err != nil {
		return Map{}, fmt.Errorf("open MANE file: %w", err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return Map{}, fmt.Errorf("gzip MANE file: %w", err)
		}
		defer gz.Close()
		r = gz
	}

	scanner := bufio.NewScanner(r)
	if !scanner.Scan() {
		return Map{}, fmt.Errorf("MANE file: empty")
	}
	header := strings.Split(scanner.Text(), "\t")
	refseqIdx := colIndex(header, "RefSeq_nuc")
	ensemblIdx := colIndex(header, "Ensembl_nuc")
	if refseqIdx < 0 || ensemblIdx < 0 {
		return Map{}, fmt.Errorf("MANE file: missing RefSeq_nuc or Ensembl_nuc column")
	}

	m := Map{
		toENST:      make(map[string]string, 20000),
		nmVersioned: make(map[string]string, 20000),
	}
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= ensemblIdx {
			continue
		}
		nmV := fields[refseqIdx]  // versioned, e.g. "NM_000546.6"
		enst := fields[ensemblIdx] // versioned, e.g. "ENST00000269305.9"
		nm := stripVersion(nmV)
		enst = stripVersion(enst)
		if nm != "" && enst != "" {
			m.toENST[nm] = enst
			m.nmVersioned[nm] = nmV
		}
	}
	return m, scanner.Err()
}

func stripVersion(s string) string {
	if dot := strings.LastIndexByte(s, '.'); dot >= 0 {
		return s[:dot]
	}
	return s
}

func colIndex(header []string, name string) int {
	for i, h := range header {
		if h == name {
			return i
		}
	}
	return -1
}
