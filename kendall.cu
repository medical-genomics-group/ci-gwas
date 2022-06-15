#include <math.h>

#include "gpuerrors.h"
#include "kendall.h"

// this will need an extern "C" declaration in the header, like the cuPC thing (this would be the
// Skeleton analogue I think) This computes correlations under the assumption that the final
// correlation matrix fits into the GPU RAM in one piece.
// This is unrealistic, but will be good for benchmarking and testing.
void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results)
{
    // This here assumes a non compressed a.
    size_t a_bytes = num_markers * num_individuals * sizeof(float);
    unsigned char *gpu_a;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_a, a_bytes));
    HANDLE_ERROR(cudaMemcpy(&gpu_a, &a, a_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    cu_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(gpu_a, num_markers, num_individuals,
                                                               gpu_results);

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_a));
    HANDLE_ERROR(cudaFree(gpu_results));
}

// TODO: this needs to be able to decode .bed binaries
// A O(n) runtime Kendall implementation for genomic marker data.
// I use unsigned char, because we will probably put bytes in here.
// I assume no NAs, x col major, blocks are comparisons for pairs for columns (e.g. x1, x2).
// What this thing needs to do foreach pair is go through it and count
// the number of occurrences of each pair.
// Then sync.
// Then compute the correlation with the counts.
__global__ void cu_marker_corr_npn(const unsigned char *a, const size_t num_markers,
                                   const size_t num_individuals, float *results)
{
    size_t tix = threadIdx.x;

    // convert linear indices into correlation matrix into (row, col) ix
    float lin_ix = blockIdx.x;
    float l = num_markers - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - lin_ix);
    size_t row = (size_t)std::floor((-b + (b * b + 4 * c)) / -2.0) + 1;
    size_t h = -(row * row) + row * (2 * l + 1);
    size_t col = lin_ix - h + row;
    size_t col_start_x = row * num_individuals;
    size_t col_start_y = col * num_individuals;

    // (0, 0) at 0
    // (0, 1) at 1
    // (0, 2) at 2
    // (1, 0) at 3
    // (1, 1) at 4
    // (1, 2) at 5
    // (2, 0) at 6
    // (2, 1) at 7
    // (2, 2) at 8
    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    // TODO: it seems stupid to jump in memory, sequential reads are probably more efficient.
    // should have ++ increment and adjust the start.
    for (size_t i = tix; i < num_individuals; i += NUMTHREADS) {
        thread_sum[(3 * a[col_start_x + i]) + a[col_start_y + i]] += 1.f;
    }

    for (size_t i = 0; i < 9; i++) {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0) {
        // produce single sum
        float sum[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++) {
            for (size_t j = 0; j < 9; j++) {
                sum[j] += thread_sums[i][j];
            }
        }
        float concordant = (sum[0] * (sum[4] + sum[5] + sum[7] + sum[8])) +
                           (sum[1] * (sum[5] + sum[8])) + (sum[3] * (sum[7] + sum[8])) +
                           (sum[4] * sum[8]);
        float discordant = (sum[1] * (sum[3] + sum[6])) +
                           (sum[2] * (sum[3] + sum[4] + sum[6] + sum[7])) + (sum[4] * sum[6]) +
                           (sum[5] * (sum[6] + sum[7]));
        float ties_x = (sum[0] * (sum[1] + sum[2])) + (sum[1] * sum[2]) +
                       (sum[3] * (sum[4] + sum[5])) + (sum[4] * sum[5]) +
                       (sum[6] * (sum[7] + sum[8]) + (sum[7] * sum[8]));
        float ties_y = (sum[0] * (sum[3] + sum[6])) + (sum[1] * (sum[4] + sum[7])) +
                       (sum[2] * (sum[5] + sum[8])) + (sum[3] * sum[6]) + (sum[4] + sum[7]) +
                       (sum[5] * sum[8]);

        size_t result_ix = col_start_x * num_markers + col_start_y;
        results[result_ix] =
            sin(M_PI / 2 * (concordant - discordant) /
                sqrt((concordant + discordant + ties_x) * (concordant + discordant + ties_y)));
    }
}