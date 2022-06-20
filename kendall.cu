#include <math.h>

#include "bed_lut.h"
#include "gpuerrors.h"
#include "kendall.h"

// Compute correlation matrix of all pairs of markers in 'a'.
// The data is not partitioned.
// 'a' is assumed to hold uncompressed genomic marker values without NaN.
void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results)
{
    // This here assumes a non compressed a.
    size_t x_bytes = num_markers * num_individuals * sizeof(unsigned char);
    unsigned char *gpu_x;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_x, x_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_x, x, x_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    cu_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(gpu_x, num_markers, num_individuals,
                                                               gpu_results);

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_x));
    HANDLE_ERROR(cudaFree(gpu_results));
}

// A O(n) runtime Kendall implementation for uncompressed genomic marker data.
__global__ void cu_marker_corr_npn(const unsigned char *x, const size_t num_markers,
                                   const size_t num_individuals, float *results)
{
    size_t tix = threadIdx.x;

    // convert linear indices into correlation matrix into (row, col) ix
    size_t lin_ix = blockIdx.x;
    float lin_ix_f = (size_t)lin_ix;
    float l = num_markers - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - lin_ix_f);
    float row = std::floor((-b + sqrt(b * b + 4 * c)) / -2.0) + 1.0;
    float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
    // offset of 1 because we don't compute the diagonal
    float col = lin_ix_f - h + row + 1;
    size_t col_start_a = row * num_individuals;
    size_t col_start_b = col * num_individuals;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];
    // TODO: it seems stupid to jump in memory, sequential reads are probably more efficient.
    // should have ++ increment and adjust the start.
    for (size_t i = tix; i < num_individuals; i += NUMTHREADS) {
        thread_sum[(3 * x[col_start_a + i]) + x[col_start_b + i]] += 1.f;
    }

    for (size_t i = 0; i < 9; i++) {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0) {
        // printf("block [x: %f; y: %f]: making single sum", row, col);

        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++) {
            for (size_t j = 0; j < 9; j++) {
                s[j] += thread_sums[i][j];
            }
        }
        float p = ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
                   (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q = ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
                   (s[5] * (s[6] + s[7])));
        float t = ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
                   (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u = ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
                   (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        // printf("linear ix: %f, row: %f, col: %f, h: %f, l: %f, corr result: %f \n", lin_ix_f,
        // row,
        //        col, h, l, kendall_corr);

        results[lin_ix] = sin(M_PI / 2 * kendall_corr);
    }
}

// Kendall correlation computation for pairs of markers compressed in .bed format
// without leading magic numbers.
// Correlation matrix is computed in one go, without any partitioning.
void cu_bed_corr_npn(const unsigned char *x, const size_t num_markers, const size_t num_individuals,
                     float *results)
{
    // this is ceil
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t x_bytes = col_len_bytes * num_markers;
    unsigned char *gpu_x;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_x, x_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_x, x, x_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    cu_bed_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(
        gpu_x, num_markers, num_individuals, col_len_bytes, gpu_results);
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_x));
    HANDLE_ERROR(cudaFree(gpu_results));
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void cu_bed_marker_corr_npn(const unsigned char *x, const size_t num_markers,
                                       const size_t num_individuals, const size_t col_len_bytes,
                                       float *results)
{
    size_t tix = threadIdx.x;

    // convert linear indices into correlation matrix into (row, col) ix
    size_t lin_ix = blockIdx.x;
    float lin_ix_f = (size_t)lin_ix;
    float l = num_markers - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - lin_ix_f);
    float row = std::floor((-b + sqrt(b * b + 4 * c)) / -2.0) + 1.0;
    float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
    // offset of 1 because we don't compute the diagonal
    float col = lin_ix_f - h + row + 1;
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    // TODO: it seems stupid to jump in memory, sequential reads are probably more efficient.
    // should have ++ increment and adjust the start.
    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS) {
        size_t aix = 4 * (size_t)(x[col_start_a + i]);
        size_t bix = 4 * (size_t)(x[col_start_b + i]);
        for (size_t j = 0; j < 4; j++) {
            if ((i * 4 + j) < num_individuals) {
                float val_a = bed_lut_a[(aix + j)];
                float val_b = bed_lut_a[(bix + j)];
                size_t comp_ix = (size_t)((3 * val_a + val_b));
                thread_sum[comp_ix] += 1.f;
            }
        }
    }

    for (size_t i = 0; i < 9; i++) {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0) {
        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++) {
            for (size_t j = 0; j < 9; j++) {
                s[j] += thread_sums[i][j];
            }
        }
        float p = ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
                   (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q = ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
                   (s[5] * (s[6] + s[7])));
        float t = ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
                   (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u = ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
                   (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        results[lin_ix] = sin(M_PI / 2 * kendall_corr);
    }
}
