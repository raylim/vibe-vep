# GPU Acceleration Feasibility Study

**Date**: March 2026  
**Hardware**: Intel Xeon Platinum 8380 (160 threads), 3× NVIDIA A100 80GB PCIe  
**CUDA**: 12.3 (nvcc), driver 580.95.05, compute cap sm_80  
**Baseline throughput**: ~720k variants/sec (CPU goroutine pool)

---

## Summary

**Codon translation and reverse complement do not benefit from GPU** at variant-annotation
batch sizes. The kernel launch overhead (~2.2 ms) exceeds the compute time for every
realistic workload.

However, **annotation source enrichment** — the sequential SQLite lookup phase that
runs after parallel annotation — is a real bottleneck at scale. Replacing it with an
in-memory open-addressing hash table yields a **656× speedup** (25,128 ns → 38 ns per
variant), and the same hash table is designed to be loaded onto the GPU for batch
processing via CUDA (see Phase 2 below).

---

## Phase 1: Codon / Revcomp GPU Spike

### Profiling Results

Benchmarks run via `go test ./internal/annotate/ -bench=. -benchmem -cpuprofile -memprofile`.

#### CPU hotspots (flat time, top 10)

| Function | Flat% | Notes |
|---|---|---|
| `runtime.memmove` | 9.4% | String copies |
| `runtime.mapaccess2_faststr` | 6.3% | **Codon table map lookup** |
| `annotate.Complement` | 5.1% | Per-base complement |
| `aeshashbody` | 3.8% | Map string hashing (codon table) |
| `runtime.slicebytetostring` | 3.2% | String allocations |
| `runtime.mallocgcTiny` | 3.4% | GC pressure |
| `annotate.GenomicToCDS` | 2.6% | CDS offset computation |
| `runtime.ReverseComplement` | 2.5% | Reverse complement |
| `annotate.formatCodonChange` | 2.2% | String formatting |
| `cache.FindExonIdx` | 1.4% | Exon linear scan |

#### Memory allocation hotspots (18.8 GB total across benchmarks)

| Function | Alloc% | Notes |
|---|---|---|
| `PredictConsequence` | **63%** | Dominant allocator — string building |
| `cdsPosRangeStr` | 11% | HGVS position string |
| `FormatHGVSp` | 9% | Protein notation string |
| `ReverseComplement` | 8% | string([]byte) conversion |
| `FormatHGVSc` | 3% | cDNA notation string |

#### Key benchmark numbers

```
BenchmarkPredictConsequence              3414 ns/op   2200 B/op   39 allocs/op
BenchmarkPredictCodingConsequence_SNV     334 ns/op    272 B/op    6 allocs/op
BenchmarkTranslateCodon                   82 ns/op      0 B/op     0 allocs/op
BenchmarkReverseComplement               127 ns/op     32 B/op     3 allocs/op
```

`PredictConsequence` processes 8 variants across exon/intron/UTR types and runs at
~293k calls/sec single-threaded. Scaled to 160 cores: ~47M calls/sec, consistent
with observed 720k variants/sec (each variant hits ~65 transcripts on average).

---

### GPU Candidate Suitability Analysis

#### 1. Codon translation

**What**: Translate 3-base codons to amino acid single-letter codes.  
**Current impl**: `map[string]byte` hash lookup — 82 ns/call, 6% of CPU time.  
**Data regularity**: High (64 possible inputs, uniform output size).  
**GPU fit**: Low — the computation per codon is 1 table lookup (trivial). Kernel
launch overhead completely dominates.

**Benchmark results** (CPU array fallback vs CUDA):

| Batch size | CPU (ns) | GPU (ns) | CPU throughput | GPU throughput |
|---|---|---|---|---|
| 1 codon | 18 | 2,180,000 | 170 MB/s | ~0 MB/s |
| 64 codons | 449 | 2,328,000 | 427 MB/s | 0.08 MB/s |
| 1k codons | 6,835 | 2,267,000 | 439 MB/s | 1.3 MB/s |
| 64k codons | 428,000 | 3,181,000 | 449 MB/s | 60 MB/s |
| **1M codons** | **6,559,000** | **2,473,000** | **457 MB/s** | **1,213 MB/s** |

GPU breaks even around **~300k codons per kernel call**. Each variant annotation
calls codon translation 1–3 times. Even batching an entire 100k-variant file yields
only ~300k codons total — right at the crossover, but with PCIe transfer overhead
making it impractical for real pipelines.

**Verdict**: ❌ Not suitable. CPU array lookup (replacing the map) is the right fix.

#### 2. Reverse complement

**What**: Complement each base and reverse the string.  
**Current impl**: Character-by-character loop, 127 ns for short sequences.  
**GPU fit**: Low for variant-length sequences (<100 bp). High only for bulk FASTA.

| Sequence length | CPU (ns) | GPU (ns) | CPU GB/s | GPU GB/s |
|---|---|---|---|---|
| 100 bp | 210 | 3,619,000 | 0.48 | 0.00003 |
| 1 kbp | 1,844 | 2,213,000 | 0.54 | 0.00045 |
| 100 kbp | 188,000 | 2,381,000 | 0.53 | 0.042 |
| **10 Mbp** | **17,536,000** | **4,637,000** | **0.57** | **2.16** |

GPU wins only at **≥3 Mbp** contiguous sequences. Variant ref/alt alleles are
typically 1–50 bp.

**Verdict**: ❌ Not suitable for current workload.

#### 3–5. Other candidates (branch logic, interval tree, pathogenicity inference)

- `PredictConsequence`: severe thread divergence → ❌
- Interval tree: pointer-chasing, fits in L3 → ❌  
- Pathogenicity inference (on-the-fly neural net): ✅ **best long-term GPU candidate**
  but requires replacing lookup-based with inference-based scoring.

---

## Phase 2: Annotation Source Enrichment (Implemented)

### Problem

`genomicindex.Lookup()` runs **sequentially** in the ordered collector callback —
one prepared-statement SQLite query per variant, typically 1–5 µs each.
At 1M variants with annotation sources enabled, this alone takes 1–5 seconds of
pure single-threaded SQLite I/O **after** the parallel annotation phase finishes.

This is the dominant source bottleneck that GPU-style batch processing directly solves.

### Solution: In-Memory Hash Table (`internal/gpuhash`)

An open-addressing hash table with 32-byte Slots, designed for GPU memory layouts:

| Field | Size | Notes |
|---|---|---|
| key hash | 8 bytes | FNV-1a of chrom+pos+ref+alt |
| am_score | 4 bytes | float32 |
| am_class | 1 byte | enum (4 values) |
| cv_sig | 1 byte | enum (25 known values) |
| cv_revstat | 1 byte | enum (9 values) |
| sig_status | 1 byte | enum (3 values) |
| sig_count | 4 bytes | uint32 |
| sig_freq | 4 bytes | float32 |
| reserved | 8 bytes | warp alignment padding |

CVClnDN (variable-length disease name) stored in a sidecar `map[uint64]string`.

**Load factor**: ≤70% → ~105M capacity slots × 32 bytes = ~3.4 GB for hg38
(AlphaMissense 72M + ClinVar 1M + SIGNAL 100k rows combined).

### Benchmark Results

```
BenchmarkLookupSQLite-160     47852    25128 ns/op   1243 B/op   40 allocs/op
BenchmarkLookupInMem-160   31129831       38 ns/op      0 B/op    0 allocs/op

Speedup: 656×  (per-variant, zero allocations)

BenchmarkBatchLookupInMem/n=1      8313457    146 ns/op
BenchmarkBatchLookupInMem/n=64      228859   5050 ns/op   (79 ns/variant)
BenchmarkBatchLookupInMem/n=1000     18106  66625 ns/op   (67 ns/variant)
BenchmarkBatchLookupInMem/n=64000      330  3.03ms ns/op   (47 ns/variant)
```

### Pipeline Integration

The `--preload-annotations` flag loads the entire genomic index into RAM at startup.
When active:

1. Per-variant `Lookup()` uses the hash table (38 ns instead of 25 µs)
2. For MAF files, `runMAFOutput` detects the preloaded source and switches to the
   **batch path**: all variants are collected, one `BatchAnnotate` call enriches them
   all at once, then rows are written — removing per-variant call overhead entirely

```bash
# Per-variant fast path (streaming, 38 ns/variant vs 25 µs)
vibe-vep annotate maf --preload-annotations input.maf

# Batch enrichment (GPU-ready: single BatchAnnotate call for entire file)
vibe-vep annotate maf --preload-annotations -o output.maf input.maf
```

**Trade-off**: ~45 s startup + ~3.4 GB RAM for hg38. Worthwhile for files >100k variants
where the source enrichment time (normally 2–5 s) exceeds the startup amortized cost.

### GPU Extension Path

The `gpuhash_cuda.cu` kernel (`internal/gpuhash/`) implements the same hash table
on GPU. When compiled with `make gpuhash CUDA_ARCH=sm_80`:

- Device-side table: 32-byte Slots, one thread per query, linear probing
- `gpuhash_batch_lookup(keys, n, results)` C export handles H2D/kernel/D2H
- `BatchLookup` on the CUDA build sends all keys in one kernel call (~1 ms launch)
- At 100k+ variants per batch, GPU lookup throughput exceeds CPU

Build commands:
```bash
make gpuhash          # compile libgpuhash.so (nvcc sm_80)
make build-gpuhash    # link vibe-vep-gpuhash with CUDA path
make bench-gpuhash    # CPU vs CUDA benchmark comparison
```

---

## Recommendations

### Short-term: CPU optimizations (implement now)

1. **Replace codon map with array lookup**: the `codonTableCPU [64]byte` in
   `internal/cudaops/cudaops_cpu.go` is the reference implementation.
   Replacing `map[string]byte` in `TranslateCodon` eliminates 10% CPU overhead.

2. **Reduce string allocations in `PredictConsequence`** (63% of all allocations):
   Use `sync.Pool` for `ConsequenceResult`, pass `[]byte` buffers.

3. **Pool `ReverseComplement` output buffers**: avoid `string(result)` heap alloc
   by returning `[]byte` where callers permit.

### Long-term: GPU inference path

If vibe-vep adds on-the-fly pathogenicity scoring (neural network), the architecture:

```
Variant stream → CPU workers (PredictConsequence) → GPU batch queue → CUDA inference → merge
```

The `internal/cudaops/` scaffold (pinned memory, CGO bridge) and `internal/gpuhash/`
(CUDA kernel build infrastructure) provide the foundation for this path.

---

## Artifacts

| Path | Description |
|---|---|
| `internal/cudaops/codon_cuda.cu` | CUDA kernels: batched codon translation + reverse complement |
| `internal/cudaops/cudaops.go` | CGO wrapper (`//go:build cuda`) |
| `internal/cudaops/cudaops_cpu.go` | CPU fallback (`//go:build !cuda`) |
| `internal/cudaops/cudaops_bench_test.go` | CPU vs GPU benchmarks across batch sizes |
| `internal/gpuhash/value.go` | 32-byte Slot layout, enum encodings |
| `internal/gpuhash/hash.go` | FNV-1a HashKey |
| `internal/gpuhash/gpuhash_cpu.go` | CPU open-addressing hash table |
| `internal/gpuhash/gpuhash_cuda.cu` | GPU batch lookup kernel |
| `internal/gpuhash/gpuhash.go` | CGO bridge (`//go:build cuda`) |
| `internal/genomicindex/inmem.go` | InMemTable, LoadInMemTable, BatchLookup, Preload |
| `internal/genomicindex/inmem_bench_test.go` | SQLite vs in-memory benchmarks |
| `internal/genomicindex/source.go` | GenomicSource.BatchAnnotate, Preload, IsPreloaded |
| `Makefile` targets `cudaops`, `gpuhash`, `build-cuda`, `build-gpuhash`, `bench-*` | Build and benchmark |


**Date**: March 2026  
**Hardware**: Intel Xeon Platinum 8380 (160 threads), 3× NVIDIA A100 80GB PCIe  
**CUDA**: 12.3 (nvcc), driver 580.95.05, compute cap sm_80  
**Baseline throughput**: ~720k variants/sec (CPU goroutine pool)

---

## Summary

**GPU acceleration does not benefit the current vibe-vep annotation workload.**

The dominant operations are branch-heavy logic and map/string operations that are
fundamentally poor fits for GPU SIMD execution. The GPU kernel launch overhead (~2.2ms)
exceeds the compute time for every realistic workload size encountered during annotation.

The best available performance wins are CPU-side and have nothing to do with the GPU:
replacing the codon map lookup with an array lookup (already implemented in `internal/cudaops`)
and reducing string allocations in `PredictConsequence`.

---

## Profiling Results

Benchmarks run via `go test ./internal/annotate/ -bench=. -benchmem -cpuprofile -memprofile`.

### CPU hotspots (flat time, top 10)

| Function | Flat% | Notes |
|---|---|---|
| `runtime.memmove` | 9.4% | String copies |
| `runtime.mapaccess2_faststr` | 6.3% | **Codon table map lookup** |
| `annotate.Complement` | 5.1% | Per-base complement |
| `aeshashbody` | 3.8% | Map string hashing (codon table) |
| `runtime.slicebytetostring` | 3.2% | String allocations |
| `runtime.mallocgcTiny` | 3.4% | GC pressure |
| `annotate.GenomicToCDS` | 2.6% | CDS offset computation |
| `runtime.ReverseComplement` | 2.5% | Reverse complement |
| `annotate.formatCodonChange` | 2.2% | String formatting |
| `cache.FindExonIdx` | 1.4% | Exon linear scan |

### Memory allocation hotspots (18.8 GB total across benchmarks)

| Function | Alloc% | Notes |
|---|---|---|
| `PredictConsequence` | **63%** | Dominant allocator — string building |
| `cdsPosRangeStr` | 11% | HGVS position string |
| `FormatHGVSp` | 9% | Protein notation string |
| `ReverseComplement` | 8% | string([]byte) conversion |
| `FormatHGVSc` | 3% | cDNA notation string |

### Key benchmark numbers

```
BenchmarkPredictConsequence              3414 ns/op   2200 B/op   39 allocs/op
BenchmarkPredictCodingConsequence_SNV     334 ns/op    272 B/op    6 allocs/op
BenchmarkTranslateCodon                   82 ns/op      0 B/op     0 allocs/op
BenchmarkReverseComplement               127 ns/op     32 B/op     3 allocs/op
```

`PredictConsequence` processes 8 variants across exon/intron/UTR types and runs at
~293k calls/sec single-threaded. Scaled to 160 cores: ~47M calls/sec, consistent
with observed 720k variants/sec (each variant hits ~65 transcripts on average).

---

## GPU Candidate Suitability Analysis

### 1. Codon translation

**What**: Translate 3-base codons to amino acid single-letter codes.  
**Current impl**: `map[string]byte` hash lookup — 82 ns/call, 6% of CPU time.  
**Data regularity**: High (64 possible inputs, uniform output size).  
**GPU fit**: Low — the computation per codon is 1 table lookup (trivial). Kernel
launch overhead completely dominates.

**Benchmark results** (CPU array fallback vs CUDA):

| Batch size | CPU (ns) | GPU (ns) | CPU throughput | GPU throughput |
|---|---|---|---|---|
| 1 codon | 18 | 2,180,000 | 170 MB/s | ~0 MB/s |
| 64 codons | 449 | 2,328,000 | 427 MB/s | 0.08 MB/s |
| 1k codons | 6,835 | 2,267,000 | 439 MB/s | 1.3 MB/s |
| 64k codons | 428,000 | 3,181,000 | 449 MB/s | 60 MB/s |
| **1M codons** | **6,559,000** | **2,473,000** | **457 MB/s** | **1,213 MB/s** |

GPU breaks even around **~300k codons per kernel call**. Each variant annotation
calls codon translation 1–3 times. Even batching an entire 100k-variant file yields
only ~300k codons total — right at the crossover, but with PCIe transfer overhead
making it impractical for real pipelines.

**Verdict**: ❌ Not suitable. CPU array lookup (replacing the map) is the right fix —
see recommendation below.

### 2. Reverse complement

**What**: Complement each base and reverse the string.  
**Current impl**: Character-by-character loop, 127 ns for short sequences.  
**GPU fit**: Low for variant-length sequences (<100 bp). High only for bulk FASTA.

| Sequence length | CPU (ns) | GPU (ns) | CPU GB/s | GPU GB/s |
|---|---|---|---|---|
| 100 bp | 210 | 3,619,000 | 0.48 | 0.00003 |
| 1 kbp | 1,844 | 2,213,000 | 0.54 | 0.00045 |
| 100 kbp | 188,000 | 2,381,000 | 0.53 | 0.042 |
| **10 Mbp** | **17,536,000** | **4,637,000** | **0.57** | **2.16** |

GPU wins only at **≥3 Mbp** contiguous sequences. Variant ref/alt alleles are
typically 1–50 bp. The CDS sequence per transcript is ~1–10 kbp — still far below
the crossover.

**Verdict**: ❌ Not suitable for current workload. GPU is useful only for whole-chromosome
FASTA operations (not in the hot path).

### 3. `PredictConsequence` (branch-heavy logic)

**What**: Determines consequence type (missense, frameshift, splice, UTR, etc.)
via a series of conditional branches across exon/intron/CDS structure.  
**GPU fit**: Very low. Consequence prediction involves deeply nested conditionals
with data-dependent paths — this causes severe **thread divergence** in SIMD GPU
execution. Divergent warps serialize execution, eliminating the throughput advantage.
The variant-specific data (exon arrays, CDS sequence, strand) lives in irregular
pointer-linked Go structs that cannot be transferred to GPU efficiently.

**Verdict**: ❌ Not suitable. Would require a complete rewrite in GPU-friendly data
layouts with no expected speedup given the branchy control flow.

### 4. Interval tree (`FindTranscripts`)

**What**: O(log n + k) augmented interval tree query per variant position.  
**GPU fit**: Low. Pointer-chasing tree traversal is the antithesis of GPU workloads.
A flat sorted-array + GPU binary search (thrust) is architecturally possible, but:
- PCIe round-trip per query adds ~10 µs vs ~0.5 µs for CPU lookup
- Batch queries possible, but transcript overlap is position-specific
- Tree already fits in L3 cache on CPU

**Verdict**: ❌ Not suitable. SIMD-friendly flat array with CPU SIMD binary search
(or even the current interval tree) outperforms GPU at this data size.

### 5. Pathogenicity inference (future — highest potential)

**What**: AlphaMissense, ClinVar lookups are currently pre-computed SQLite
point-lookups (~1–5 µs each, sequential). If replaced with **batched neural
network inference** (e.g., a re-implemented AlphaMissense model), the GPU would
be highly effective.  
**GPU fit**: High — small neural networks run at tens of millions of samples/sec
on an A100.

**Verdict**: ✅ **Best GPU candidate** — but requires workload change from lookup
to inference. Out of scope for current architecture.

---

## Recommendations

### Short-term: CPU optimizations (implement now)

These are straightforward changes derived from the profiling data:

1. **Replace codon map with array lookup** (`internal/annotate/codon.go`):
   The `codonTable map[string]byte` causes 10% of CPU time via map hashing
   (`mapaccess2_faststr` + `aeshashbody`). The `internal/cudaops` package already
   implements the correct array-indexed lookup (`codonTableCPU [64]byte`).
   Replacing the map in `TranslateCodon` with 2-bit base encoding + array indexing
   would eliminate that 10% overhead with no correctness risk.

2. **Reduce string allocations in `PredictConsequence`** (63% of all allocations):
   Use `sync.Pool` for `ConsequenceResult`, pass `[]byte` buffers instead of
   building `string` values for codon changes, or use a stack-allocated arena for
   the per-variant result. This would reduce GC pressure significantly.

3. **Pool `ReverseComplement` output buffers**: The function already uses a
   stack-allocated `[64]byte` for the complement but still calls `string(result)`
   which heap-allocates. Consider changing the return type to `[]byte` where
   callers allow it.

### Long-term: GPU integration path (only if workload evolves)

If vibe-vep evolves to include on-the-fly pathogenicity scoring (rather than
pre-computed lookup), the following architecture would work:

```
Variant stream → CPU workers (PredictConsequence) → GPU batch queue → CUDA inference kernel → merge results
```

- Use a dedicated GPU worker goroutine with a batch accumulator
- Fill a pinned-memory buffer until batch size ≥ 1k variants or a timeout fires
- Invoke CUDA inference kernel asynchronously
- Merge GPU scores with CPU annotation results in the ordered collector
- See `internal/cudaops/` for the CGO/pinned-memory infrastructure scaffold

The `internal/cudaops/` package provides the build infrastructure, CGO bridge,
pinned-memory allocator, and benchmarks needed to extend this path.

---

## Artifacts

| Path | Description |
|---|---|
| `internal/cudaops/codon_cuda.cu` | CUDA kernels: batched codon translation + reverse complement |
| `internal/cudaops/cudaops.go` | Go CGO wrapper (`//go:build cuda`) |
| `internal/cudaops/cudaops_cpu.go` | CPU fallback (`//go:build !cuda`) |
| `internal/cudaops/cudaops_bench_test.go` | CPU vs GPU benchmarks across batch sizes |
| `Makefile` targets `cudaops`, `build-cuda`, `bench-cuda` | Build and benchmark commands |

### Build commands

```bash
# Build the CUDA shared library (requires nvcc)
make cudaops CUDA_ARCH=sm_80

# Run CPU vs GPU comparison benchmarks
make bench-cuda

# Build main binary with CUDA enabled
make build-cuda
```
