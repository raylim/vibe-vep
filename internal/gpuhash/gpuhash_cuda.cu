/*
 * gpuhash_cuda.cu – GPU-resident open-addressing hash table for batch
 * genomic annotation lookups (AlphaMissense, ClinVar, SIGNAL).
 *
 * Layout
 * ------
 * Each slot is 32 bytes (natural alignment for 32-bit warp lanes):
 *
 *   offset  size  field
 *    0       8    key_hash (uint64_t; 0 = empty)
 *    8       4    am_score (float)
 *   12       1    am_class (uint8_t)
 *   13       1    cv_sig   (uint8_t)
 *   14       1    cv_revstat (uint8_t)
 *   15       1    sig_status (uint8_t)
 *   16       4    sig_count  (uint32_t)
 *   20       4    sig_freq   (float)
 *   24       8    reserved
 *
 * Kernel
 * ------
 * batch_lookup_kernel: one thread per query.  Linear probing; each thread
 * walks its probe chain independently with no inter-thread communication.
 * Result slots map 1-to-1 with the input hash array: not-found → zero-filled.
 *
 * Concurrency
 * -----------
 * The table is built once on the host and copied to device memory before any
 * lookups.  After the copy the device-side table is read-only, so concurrent
 * kernel invocations are safe.
 */

#include <stdint.h>
#include <string.h>

/* ── slot layout ──────────────────────────────────────────────── */

typedef struct __attribute__((packed, aligned(8))) {
    uint64_t hash;
    float    am_score;
    uint8_t  am_class;
    uint8_t  cv_sig;
    uint8_t  cv_revstat;
    uint8_t  sig_status;
    uint32_t sig_count;
    float    sig_freq;
    uint8_t  reserved[8];
} Slot; /* 32 bytes */

/* ── result layout (mirrors Go Value fields the GPU fills) ──── */

typedef struct __attribute__((packed, aligned(4))) {
    float    am_score;
    uint8_t  am_class;
    uint8_t  cv_sig;
    uint8_t  cv_revstat;
    uint8_t  sig_status;
    uint32_t sig_count;
    float    sig_freq;
    uint8_t  found;   /* 1 = hit, 0 = miss */
    uint8_t  pad[3];
} Result; /* 20 bytes */

/* ── CUDA kernel ─────────────────────────────────────────────── */

__global__ void batch_lookup_kernel(
    const Slot   * __restrict__ table,
    uint64_t                    capacity,   /* power of two */
    const uint64_t * __restrict__ hashes,
    Result       * __restrict__ results,
    uint32_t                    n_queries)
{
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_queries) return;

    uint64_t h = hashes[tid];
    Result   r;
    memset(&r, 0, sizeof(r));

    if (h != 0) {
        uint64_t mask = capacity - 1;
        uint64_t idx  = h & mask;
        for (uint64_t probe = 0; probe < capacity; ++probe) {
            const Slot *s = table + idx;
            if (s->hash == 0) break;          /* empty slot → miss */
            if (s->hash == h) {               /* hit */
                r.am_score   = s->am_score;
                r.am_class   = s->am_class;
                r.cv_sig     = s->cv_sig;
                r.cv_revstat = s->cv_revstat;
                r.sig_status = s->sig_status;
                r.sig_count  = s->sig_count;
                r.sig_freq   = s->sig_freq;
                r.found      = 1;
                break;
            }
            idx = (idx + 1) & mask;
        }
    }

    results[tid] = r;
}

/* ── C-callable exports ─────────────────────────────────────── */

extern "C" {

/*
 * gpuhash_create – allocate a device-side hash table.
 * capacity must be a power of two.
 * Returns pointer to device memory (cast to void*), or NULL on error.
 */
void *gpuhash_create(uint64_t capacity) {
    Slot *d_table = nullptr;
    size_t bytes = (size_t)capacity * sizeof(Slot);
    if (cudaMalloc(&d_table, bytes) != cudaSuccess) return nullptr;
    cudaMemset(d_table, 0, bytes);
    return d_table;
}

/*
 * gpuhash_load – copy host slots into the device table.
 * d_table: device pointer returned by gpuhash_create.
 * h_slots: host array of Slot structs.
 * count:   number of slots to copy.
 * offset:  first slot index to write (for chunked loading).
 */
int gpuhash_load(void *d_table, const Slot *h_slots, uint64_t count, uint64_t offset) {
    Slot *dst = (Slot *)d_table + offset;
    if (cudaMemcpy(dst, h_slots, count * sizeof(Slot), cudaMemcpyHostToDevice) != cudaSuccess)
        return -1;
    return 0;
}

/*
 * gpuhash_destroy – free device memory.
 */
void gpuhash_destroy(void *d_table) {
    if (d_table) cudaFree(d_table);
}

/*
 * gpuhash_batch_lookup – launch batch lookup kernel.
 * d_table:    device hash table pointer.
 * capacity:   number of slots (power of two).
 * d_hashes:   device array of n uint64_t query hashes.
 * d_results:  device array of n Result structs (output).
 * n:          number of queries.
 * Returns 0 on success, -1 on CUDA error.
 */
int gpuhash_batch_lookup(
    void     *d_table,
    uint64_t  capacity,
    void     *d_hashes,
    void     *d_results,
    uint32_t  n)
{
    if (n == 0) return 0;
    int block = 256;
    int grid  = (n + block - 1) / block;
    batch_lookup_kernel<<<grid, block>>>(
        (const Slot *)d_table,
        capacity,
        (const uint64_t *)d_hashes,
        (Result *)d_results,
        n);
    return (cudaGetLastError() == cudaSuccess) ? 0 : -1;
}

/*
 * gpuhash_sync – wait for all device operations to complete.
 */
int gpuhash_sync(void) {
    return (cudaDeviceSynchronize() == cudaSuccess) ? 0 : -1;
}

/*
 * gpuhash_alloc_device – allocate raw device memory (for hashes/results buffers).
 */
void *gpuhash_alloc_device(size_t bytes) {
    void *p = nullptr;
    return (cudaMalloc(&p, bytes) == cudaSuccess) ? p : nullptr;
}

/*
 * gpuhash_free_device – free device memory allocated by gpuhash_alloc_device.
 */
void gpuhash_free_device(void *p) {
    if (p) cudaFree(p);
}

/*
 * gpuhash_copy_hashes_to_device – copy host hash array to device buffer.
 */
int gpuhash_copy_hashes_to_device(void *d_buf, const uint64_t *h_hashes, uint32_t n) {
    return (cudaMemcpy(d_buf, h_hashes, (size_t)n * sizeof(uint64_t), cudaMemcpyHostToDevice)
                == cudaSuccess) ? 0 : -1;
}

/*
 * gpuhash_copy_results_to_host – copy result array from device to host.
 * result_bytes: total bytes = n * sizeof(Result) = n * 20
 */
int gpuhash_copy_results_to_host(void *h_buf, const void *d_buf, size_t result_bytes) {
    return (cudaMemcpy(h_buf, d_buf, result_bytes, cudaMemcpyDeviceToHost)
                == cudaSuccess) ? 0 : -1;
}

/*
 * gpuhash_result_size – returns sizeof(Result) so CGO can allocate correctly.
 */
size_t gpuhash_result_size(void) { return sizeof(Result); }

} /* extern "C" */
