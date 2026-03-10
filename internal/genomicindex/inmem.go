package genomicindex

import (
	"fmt"
	"math"
	"strconv"

	"github.com/inodb/vibe-vep/internal/gpuhash"
)

// InMemTable is a gpuhash.Table loaded from the genomic_annotations SQLite
// database.  It provides O(1) hash-table lookups instead of prepared-statement
// executions, removing the sequential bottleneck in the annotation pipeline.
//
// The table is read-only once built; all methods are safe for concurrent use.
type InMemTable struct {
	tbl *gpuhash.Table
	// clnDN maps hash → disease name string for ClinVar entries.
	// We store it separately because it is variable-length and rarely needed.
	clnDN map[uint64]string
}

// LoadInMemTable reads every row from the genomic_annotations table in s and
// builds a gpuhash.Table plus a sidecar clnDN map.
//
// For a typical hg38 database (AlphaMissense ~72 M rows, ClinVar ~1 M rows,
// SIGNAL ~100 k rows) this takes ~45 seconds and uses ~3 GB of RAM.  The
// caller should keep the Store alive for the process lifetime.
//
// After LoadInMemTable returns without error, single-variant Lookup calls are
// served from the in-memory table transparently.
func (s *Store) LoadInMemTable() (*InMemTable, error) {
	// Count rows to pre-size the table.
	var rowCount int64
	err := s.db.QueryRow("SELECT COUNT(*) FROM genomic_annotations").Scan(&rowCount)
	if err != nil {
		return nil, fmt.Errorf("count rows: %w", err)
	}

	// Allocate with 43% headroom so load factor ≤ 0.70.
	capacity := uint64(math.Ceil(float64(rowCount) / 0.70))
	tbl, err := gpuhash.NewTable(capacity)
	if err != nil {
		return nil, fmt.Errorf("new table: %w", err)
	}

	clnDN := make(map[uint64]string, rowCount/100) // most ClinVar entries have a disease name

	rows, err := s.db.Query(`SELECT chrom, pos, ref, alt,
		am_score, am_class,
		cv_clnsig, cv_revstat, cv_clndn,
		sig_mut_status, sig_count, sig_freq
		FROM genomic_annotations`)
	if err != nil {
		tbl.Free()
		return nil, fmt.Errorf("query all rows: %w", err)
	}
	defer rows.Close()

	var (
		chrom, ref, alt                           string
		amClass, cvSig, cvRevStat, cvClnDN        string
		sigMut, sigCount, sigFreq                 string
		pos                                       int64
		amScore                                   float64
	)

	for rows.Next() {
		if err := rows.Scan(
			&chrom, &pos, &ref, &alt,
			&amScore, &amClass,
			&cvSig, &cvRevStat, &cvClnDN,
			&sigMut, &sigCount, &sigFreq,
		); err != nil {
			tbl.Free()
			return nil, fmt.Errorf("scan row: %w", err)
		}

		h := gpuhash.HashKey(chrom, pos, ref, alt)

		sc := uint32(0)
		if sigCount != "" {
			n, _ := strconv.ParseInt(sigCount, 10, 32)
			sc = uint32(n)
		}
		sf := float32(0)
		if sigFreq != "" {
			f, _ := strconv.ParseFloat(sigFreq, 32)
			sf = float32(f)
		}

		tbl.Insert(gpuhash.Slot{
			Hash:      h,
			AMScore:   float32(amScore),
			AMClass:   gpuhash.EncodeAMClass(amClass),
			CVSig:     gpuhash.EncodeCVSig(cvSig),
			CVRevStat: gpuhash.EncodeCVRevStat(cvRevStat),
			SigStatus: gpuhash.EncodeSigStatus(sigMut),
			SigCount:  sc,
			SigFreq:   sf,
		})

		if cvClnDN != "" {
			clnDN[h] = cvClnDN
		}
	}
	if err := rows.Err(); err != nil {
		tbl.Free()
		return nil, fmt.Errorf("scan rows: %w", err)
	}

	return &InMemTable{tbl: tbl, clnDN: clnDN}, nil
}

// Preload calls LoadInMemTable and attaches the result to this Store so
// subsequent Lookup calls use the in-memory hash table automatically.
// It is safe to call Preload once at startup; concurrent Lookup calls must
// not overlap with Preload.
func (s *Store) Preload() error {
	m, err := s.LoadInMemTable()
	if err != nil {
		return err
	}
	s.inmem = m
	return nil
}

// LookupKey is a single annotation lookup request.
type LookupKey struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
}

// BatchLookup resolves all keys at once using the in-memory hash table and
// returns a Result slice in the same order.  The CVClnDN field is populated
// from the sidecar map.  Not-found entries are zero-valued Results.
func (m *InMemTable) BatchLookup(keys []LookupKey) []Result {
	if len(keys) == 0 {
		return nil
	}

	hashes := make([]uint64, len(keys))
	for i, k := range keys {
		hashes[i] = gpuhash.HashKey(k.Chrom, k.Pos, k.Ref, k.Alt)
	}

	values := make([]gpuhash.Value, len(keys))
	m.tbl.BatchLookup(hashes, values)

	results := make([]Result, len(keys))
	for i, v := range values {
		if !v.HasAM() && !v.HasCV() && !v.HasSig() {
			continue
		}
		r := &results[i]
		if v.HasAM() {
			r.AMScore = v.AMScore
			r.AMClass = gpuhash.AMClassString[v.AMClass]
		}
		if v.HasCV() {
			r.CVClnSig = gpuhash.CVSigStrings[v.CVSig]
			r.CVClnRevStat = gpuhash.CVRevStatStrings[v.CVRevStat]
			if dn, ok := m.clnDN[hashes[i]]; ok {
				r.CVClnDN = dn
			}
		}
		if v.HasSig() {
			r.SigMutStatus = gpuhash.SigStatusStrings[v.SigStatus]
			if v.SigCount > 0 {
				r.SigCount = strconv.Itoa(int(v.SigCount))
			}
			if v.SigFreq > 0 {
				r.SigFreq = strconv.FormatFloat(float64(v.SigFreq), 'f', 6, 32)
			}
		}
	}
	return results
}

// Free releases the memory held by this InMemTable.
func (m *InMemTable) Free() {
	if m.tbl != nil {
		m.tbl.Free()
	}
}

// Len returns the number of entries in the hash table.
func (m *InMemTable) Len() uint64 {
	if m.tbl == nil {
		return 0
	}
	return m.tbl.Len()
}

// lookupInMemResult performs a single lookup against the in-memory table and
// returns a Result.  Used as a faster drop-in for Store.Lookup when the
// in-memory table is available.
func (m *InMemTable) lookupInMemResult(chrom string, pos int64, ref, alt string) (Result, bool) {
	h := gpuhash.HashKey(chrom, pos, ref, alt)
	v, ok := m.tbl.Lookup(h)
	if !ok {
		return Result{}, false
	}
	r := Result{}
	if v.HasAM() {
		r.AMScore = v.AMScore
		r.AMClass = gpuhash.AMClassString[v.AMClass]
	}
	if v.HasCV() {
		r.CVClnSig = gpuhash.CVSigStrings[v.CVSig]
		r.CVClnRevStat = gpuhash.CVRevStatStrings[v.CVRevStat]
		if dn, ok := m.clnDN[h]; ok {
			r.CVClnDN = dn
		}
	}
	if v.HasSig() {
		r.SigMutStatus = gpuhash.SigStatusStrings[v.SigStatus]
		if v.SigCount > 0 {
			r.SigCount = strconv.Itoa(int(v.SigCount))
		}
		if v.SigFreq > 0 {
			r.SigFreq = strconv.FormatFloat(float64(v.SigFreq), 'f', 6, 32)
		}
	}
	return r, true
}
