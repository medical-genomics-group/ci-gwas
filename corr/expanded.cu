#include <math.h>

#include "../cuPC/gpuerrors.h"
#include "bed_lut.h"
#include "expanded.h"

// Compute correlation matrix of all pairs of markers in 'x'.
// The data is not partitioned.
// 'x' is assumed to hold uncompressed genomic marker values without NaN.
void cu_corr_npn(const unsigned char *marker_vals, const size_t num_markers,
                 const size_t num_individuals, float *results)
{
    // This here assumes a non compressed a.
    size_t marker_vals_bytes = num_markers * num_individuals * sizeof(unsigned char);
    unsigned char *gpu_marker_vals;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, marker_vals_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    cu_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(gpu_marker_vals, num_markers,
                                                               num_individuals, gpu_results);

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_results));
}

// A O(n) runtime Kendall implementation for uncompressed genomic marker data.
__global__ void cu_marker_corr_npn(const unsigned char *marker_vals, const size_t num_markers,
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
        thread_sum[(3 * marker_vals[col_start_a + i]) + marker_vals[col_start_b + i]] += 1.f;
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