package parquet

import (
	"fmt"
	"io"
	"sort"

	pq "github.com/parquet-go/parquet-go"
)

// DefaultRowGroupSize is the default number of rows per row group.
// ~20K rows gives ~20-50KB per group for typical annotation data,
// ideal for DuckDB-WASM HTTP range requests.
const DefaultRowGroupSize = 20_000

// Writer writes variant annotation rows to a Parquet file.
type Writer struct {
	w            *pq.GenericWriter[Row]
	rowGroupSize int
}

// NewWriter creates a Parquet writer with Zstd compression.
// If rowGroupSize is 0, DefaultRowGroupSize is used.
func NewWriter(w io.Writer, rowGroupSize int) *Writer {
	if rowGroupSize <= 0 {
		rowGroupSize = DefaultRowGroupSize
	}
	schema := pq.SchemaOf(Row{})
	pw := pq.NewGenericWriter[Row](w,
		schema,
		pq.Compression(&pq.Zstd),
		pq.CreatedBy("vibe-vep", "", ""),
	)
	return &Writer{w: pw, rowGroupSize: rowGroupSize}
}

// WriteRows writes a batch of rows. Rows should be pre-sorted.
func (w *Writer) WriteRows(rows []Row) error {
	for i := 0; i < len(rows); i += w.rowGroupSize {
		end := i + w.rowGroupSize
		if end > len(rows) {
			end = len(rows)
		}
		if _, err := w.w.Write(rows[i:end]); err != nil {
			return fmt.Errorf("writing row group: %w", err)
		}
	}
	return nil
}

// Close flushes and closes the Parquet writer.
func (w *Writer) Close() error {
	return w.w.Close()
}

// SortRows sorts rows by (chrom_numeric, pos, ref, alt, transcript_id)
// for optimal row group stats pruning.
func SortRows(rows []Row) {
	sort.Slice(rows, func(i, j int) bool {
		a, b := rows[i], rows[j]
		if a.ChromNumeric != b.ChromNumeric {
			return a.ChromNumeric < b.ChromNumeric
		}
		if a.Pos != b.Pos {
			return a.Pos < b.Pos
		}
		if a.Ref != b.Ref {
			return a.Ref < b.Ref
		}
		if a.Alt != b.Alt {
			return a.Alt < b.Alt
		}
		return a.TranscriptID < b.TranscriptID
	})
}
