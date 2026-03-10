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

// Map maps an unversioned RefSeq accession (e.g. "NM_000546") to its
// unversioned Ensembl transcript (e.g. "ENST00000269305").
type Map map[string]string

// ENST returns the unversioned Ensembl transcript for a RefSeq accession.
// The refseq argument may be versioned (e.g. "NM_000546.6") or not.
// Returns ("", false) if not present.
func (m Map) ENST(refseq string) (string, bool) {
	k := stripVersion(refseq)
	v, ok := m[k]
	return v, ok
}

// HasRefSeq reports whether the RefSeq accession (versioned or not) is in MANE.
func (m Map) HasRefSeq(refseq string) bool {
	_, ok := m[stripVersion(refseq)]
	return ok
}

// Load reads a MANE summary file (plain or gzipped) and returns the mapping.
// The file is the NCBI MANE summary (e.g. MANE.GRCh38.v1.5.summary.txt.gz).
func Load(path string) (Map, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open MANE file: %w", err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("gzip MANE file: %w", err)
		}
		defer gz.Close()
		r = gz
	}

	scanner := bufio.NewScanner(r)
	if !scanner.Scan() {
		return nil, fmt.Errorf("MANE file: empty")
	}
	header := strings.Split(scanner.Text(), "\t")
	refseqIdx := colIndex(header, "RefSeq_nuc")
	ensemblIdx := colIndex(header, "Ensembl_nuc")
	if refseqIdx < 0 || ensemblIdx < 0 {
		return nil, fmt.Errorf("MANE file: missing RefSeq_nuc or Ensembl_nuc column")
	}

	m := make(Map, 20000)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= ensemblIdx {
			continue
		}
		nm := stripVersion(fields[refseqIdx])
		enst := stripVersion(fields[ensemblIdx])
		if nm != "" && enst != "" {
			m[nm] = enst
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
