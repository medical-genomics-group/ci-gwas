#include <math.h>

#define NUMTHREADS 256

void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results)
{
    size_t a_bytes = num_markers * num_individuals * sizeof(float);
    float *gpu_a;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    dim3 grid(num_markers, num_markers);

    cudaMalloc(&gpu_a, a_bytes);
    // TODO: this is not going to work if a is larger than the device memory.
    // If that is the case, we need to split up a into blocks in some clever way.
    // This is what Mahdi said he solved in his R code, will look at that.
    cudaMemcpy(gpu_a, a, a_bytes, cudaMemcpyHostToDevice);
}

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
    size_t i;
    size_t j;
    size_t tests;
    size_t tix = threadIdx.x;
    size_t col_ix_x = blockIdx.x;
    size_t col_ix_y = blockIdx.y;
    size_t col_start_x = col_ix_x * num_individuals;
    size_t col_start_y = col_ix_y * num_individuals;

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
    for (i = tix; i < num_individuals; i += NUMTHREADS) {
        thread_sum[(3 * a[col_start_x + i]) + a[col_start_y + i]] += 1.f;
    }

    thread_sums[tix] = thread_sum;

    // consolidate thread_sums
    __syncthreads();
    if (tx == 0) {
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