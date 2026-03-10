// CUDA kernels for batched genomic sequence operations.
//
// Build with:
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC -shared -o libcudaops.so codon_cuda.cu
#include <cuda_runtime.h>
#include <stdint.h>

// base2idx converts a DNA base character to a 2-bit index.
// A=0, C=1, G=2, T=3. Unknown bases map to 0.
__device__ __forceinline__ uint8_t base2idx(char b) {
    switch (b) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default:             return 0;
    }
}

// codon_index encodes a 3-base codon as a 6-bit index (0–63).
__device__ __forceinline__ uint8_t codon_index(char b0, char b1, char b2) {
    return (base2idx(b0) << 4) | (base2idx(b1) << 2) | base2idx(b2);
}

// Standard genetic code indexed by codon_index(b0,b1,b2).
// Entries are single-letter amino acid codes; '*' for stop, 'X' for unknown.
//
// Index formula: (base2idx(b0) << 4) | (base2idx(b1) << 2) | base2idx(b2)
// with A=0, C=1, G=2, T=3.
//
// Row by first base:
//   A (0x00–0x0F): K N K N  T T T T  R S R S  I I M I
//   C (0x10–0x1F): Q H Q H  P P P P  R R R R  L L L L
//   G (0x20–0x2F): E D E D  A A A A  G G G G  V V V V
//   T (0x30–0x3F): * Y * Y  S S S S  * C W C  L F L F
__constant__ char kCodonTable[64] = {
    'K','N','K','N', 'T','T','T','T', 'R','S','R','S', 'I','I','M','I',
    'Q','H','Q','H', 'P','P','P','P', 'R','R','R','R', 'L','L','L','L',
    'E','D','E','D', 'A','A','A','A', 'G','G','G','G', 'V','V','V','V',
    '*','Y','*','Y', 'S','S','S','S', '*','C','W','C', 'L','F','L','F',
};

// kComplementTable maps each ASCII byte to its complement.
// Only A,C,G,T (upper/lower) are mapped; everything else maps to 'N'.
__constant__ char kComplementTable[256];

// cuda_translate_codons translates n codons in parallel.
// in_codons: packed input, 3 bytes per codon, length = 3*n.
// out_aa:    output amino acid array, length = n.
__global__ void cuda_translate_codons(
    const char* __restrict__ in_codons,
    char*       __restrict__ out_aa,
    int                      n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    const char* c = in_codons + idx * 3;
    out_aa[idx] = kCodonTable[codon_index(c[0], c[1], c[2])];
}

// cuda_reverse_complement reverse-complements a DNA sequence in parallel.
// Each thread handles one base of the output.
// out[i] = complement(in[n-1-i])
__global__ void cuda_reverse_complement(
    const char* __restrict__ in,
    char*       __restrict__ out,
    int                      n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    out[idx] = kComplementTable[(unsigned char)in[n - 1 - idx]];
}

// ---- Host-callable C functions for CGO ----

extern "C" {

// InitComplementTable fills the constant complement table.
// Must be called once before using cuda_reverse_complement.
void cudaops_init_complement_table() {
    char h[256];
    for (int i = 0; i < 256; i++) h[i] = 'N';
    h['A'] = 'T'; h['T'] = 'A'; h['G'] = 'C'; h['C'] = 'G';
    h['a'] = 't'; h['t'] = 'a'; h['g'] = 'c'; h['c'] = 'g';
    cudaMemcpyToSymbol(kComplementTable, h, 256);
}

// cudaops_translate_codons translates n codons.
// in_codons: host pointer, 3*n bytes. out_aa: host pointer, n bytes.
// Returns 0 on success, non-zero on CUDA error.
int cudaops_translate_codons(const char* in_codons, char* out_aa, int n) {
    if (n <= 0) return 0;
    char *d_in = nullptr, *d_out = nullptr;
    cudaError_t err;
    if ((err = cudaMalloc(&d_in, 3 * n)) != cudaSuccess) return (int)err;
    if ((err = cudaMalloc(&d_out, n)) != cudaSuccess) { cudaFree(d_in); return (int)err; }
    if ((err = cudaMemcpy(d_in, in_codons, 3 * n, cudaMemcpyHostToDevice)) != cudaSuccess)
        goto fail;
    {
        int threads = 256;
        int blocks = (n + threads - 1) / threads;
        cuda_translate_codons<<<blocks, threads>>>(d_in, d_out, n);
        if ((err = cudaGetLastError()) != cudaSuccess) goto fail;
        if ((err = cudaDeviceSynchronize()) != cudaSuccess) goto fail;
    }
    err = cudaMemcpy(out_aa, d_out, n, cudaMemcpyDeviceToHost);
fail:
    cudaFree(d_in);
    cudaFree(d_out);
    return (int)err;
}

// cudaops_translate_codons_pinned translates n codons using pinned host memory
// for lower PCIe latency. Caller must provide CUDA-pinned buffers.
int cudaops_translate_codons_pinned(const char* in_codons, char* out_aa, int n) {
    if (n <= 0) return 0;
    char *d_in = nullptr, *d_out = nullptr;
    cudaError_t err;
    if ((err = cudaMalloc(&d_in, 3 * n)) != cudaSuccess) return (int)err;
    if ((err = cudaMalloc(&d_out, n)) != cudaSuccess) { cudaFree(d_in); return (int)err; }
    // Async copy from pinned memory
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    if ((err = cudaMemcpyAsync(d_in, in_codons, 3*n, cudaMemcpyHostToDevice, stream)) != cudaSuccess)
        goto fail_stream;
    {
        int threads = 256;
        int blocks = (n + threads - 1) / threads;
        cuda_translate_codons<<<blocks, threads, 0, stream>>>(d_in, d_out, n);
        if ((err = cudaGetLastError()) != cudaSuccess) goto fail_stream;
        if ((err = cudaMemcpyAsync(out_aa, d_out, n, cudaMemcpyDeviceToHost, stream)) != cudaSuccess)
            goto fail_stream;
        err = cudaStreamSynchronize(stream);
    }
fail_stream:
    cudaStreamDestroy(stream);
    cudaFree(d_in);
    cudaFree(d_out);
    return (int)err;
}

// cudaops_reverse_complement reverse-complements a sequence.
// in/out are host pointers of length n.
int cudaops_reverse_complement(const char* in, char* out, int n) {
    if (n <= 0) return 0;
    char *d_in = nullptr, *d_out = nullptr;
    cudaError_t err;
    if ((err = cudaMalloc(&d_in, n)) != cudaSuccess) return (int)err;
    if ((err = cudaMalloc(&d_out, n)) != cudaSuccess) { cudaFree(d_in); return (int)err; }
    if ((err = cudaMemcpy(d_in, in, n, cudaMemcpyHostToDevice)) != cudaSuccess) goto fail;
    {
        int threads = 256;
        int blocks = (n + threads - 1) / threads;
        cuda_reverse_complement<<<blocks, threads>>>(d_in, d_out, n);
        if ((err = cudaGetLastError()) != cudaSuccess) goto fail;
        if ((err = cudaDeviceSynchronize()) != cudaSuccess) goto fail;
    }
    err = cudaMemcpy(out, d_out, n, cudaMemcpyDeviceToHost);
fail:
    cudaFree(d_in);
    cudaFree(d_out);
    return (int)err;
}

// cudaops_alloc_pinned allocates n bytes of pinned host memory.
// Returns pointer on success, NULL on failure.
void* cudaops_alloc_pinned(int n) {
    void* ptr = nullptr;
    cudaMallocHost(&ptr, n);
    return ptr;
}

// cudaops_free_pinned frees pinned host memory.
void cudaops_free_pinned(void* ptr) {
    cudaFreeHost(ptr);
}

} // extern "C"
