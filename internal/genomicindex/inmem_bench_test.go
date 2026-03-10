package genomicindex

import (
	"database/sql"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"testing"

	_ "modernc.org/sqlite"
)

// setupBenchDB creates a small SQLite database and returns an open Store.
func setupBenchDB(b *testing.B) *Store {
	b.Helper()
	dir, err := os.MkdirTemp("", "genomicindex-bench-*")
	if err != nil {
		b.Fatal(err)
	}
	b.Cleanup(func() { os.RemoveAll(dir) })

	dbPath := filepath.Join(dir, "bench.sqlite")
	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		b.Fatal(err)
	}
	_, err = db.Exec(`CREATE TABLE genomic_annotations (
		chrom TEXT NOT NULL,
		pos INTEGER NOT NULL,
		ref TEXT NOT NULL,
		alt TEXT NOT NULL,
		am_score REAL NOT NULL DEFAULT 0,
		am_class TEXT NOT NULL DEFAULT '',
		cv_clnsig TEXT NOT NULL DEFAULT '',
		cv_revstat TEXT NOT NULL DEFAULT '',
		cv_clndn TEXT NOT NULL DEFAULT '',
		sig_mut_status TEXT NOT NULL DEFAULT '',
		sig_count TEXT NOT NULL DEFAULT '',
		sig_freq TEXT NOT NULL DEFAULT '',
		PRIMARY KEY (chrom, pos, ref, alt)
	) WITHOUT ROWID`)
	if err != nil {
		b.Fatal(err)
	}
	// Insert a handful of real-ish rows.
	rows := []struct{ chrom, ref, alt string; pos int64; score float64; class string }{
		{"12", "C", "A", 25245350, 0.9876, "likely_pathogenic"},
		{"7", "G", "T", 55191822, 0.2345, "benign"},
		{"17", "A", "G", 7674220, 0.7654, "pathogenic"},
	}
	for _, r := range rows {
		_, err = db.Exec(`INSERT INTO genomic_annotations (chrom,pos,ref,alt,am_score,am_class) VALUES (?,?,?,?,?,?)`,
			r.chrom, r.pos, r.ref, r.alt, r.score, r.class)
		if err != nil {
			b.Fatal(err)
		}
	}
	db.Close()

	store, err := Open(dbPath)
	if err != nil {
		b.Fatal(err)
	}
	b.Cleanup(func() { store.Close() })
	return store
}

// BenchmarkLookupSQLite measures per-variant SQLite prepared-statement lookup.
func BenchmarkLookupSQLite(b *testing.B) {
	store := setupBenchDB(b)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = store.Lookup("1", int64(i%1_000_000), "A", "T")
	}
}

// BenchmarkLookupInMem measures per-variant in-memory hash-table lookup after
// a Preload call.
func BenchmarkLookupInMem(b *testing.B) {
	store := setupBenchDB(b)
	if err := store.Preload(); err != nil {
		b.Fatal(err)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = store.Lookup("1", int64(i%1_000_000), "A", "T")
	}
}

// BenchmarkBatchLookupInMem measures InMemTable.BatchLookup across varying batch sizes.
func BenchmarkBatchLookupInMem(b *testing.B) {
	store := setupBenchDB(b)
	if err := store.Preload(); err != nil {
		b.Fatal(err)
	}

	chroms := []string{"1", "2", "7", "12", "17"}
	rng := rand.New(rand.NewSource(42))

	for _, sz := range []int{1, 64, 1_000, 64_000} {
		sz := sz
		keys := make([]LookupKey, sz)
		for i := range keys {
			keys[i] = LookupKey{
				Chrom: chroms[rng.Intn(len(chroms))],
				Pos:   int64(rng.Intn(250_000_000)),
				Ref:   "A",
				Alt:   "T",
			}
		}

		b.Run(fmt.Sprintf("n=%d", sz), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = store.inmem.BatchLookup(keys)
			}
		})
	}
}
