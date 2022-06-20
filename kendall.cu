#include <math.h>

#include "bed_lut.h"
#include "gpuerrors.h"
#include "kendall.h"

// This computes correlations under the assumption that the final
// correlation matrix fits into the GPU RAM in one piece.
// This is unrealistic, but will be good for benchmarking and testing.
void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results)
{
    // This here assumes a non compressed a.
    size_t a_bytes = num_markers * num_individuals * sizeof(unsigned char);
    unsigned char *gpu_a;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_a, a_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_a, a, a_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    cu_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(gpu_a, num_markers, num_individuals,
                                                               gpu_results);

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_a));
    HANDLE_ERROR(cudaFree(gpu_results));
}

// A O(n) runtime Kendall implementation for uncompressed genomic marker data.
__global__ void cu_marker_corr_npn(const unsigned char *a, const size_t num_markers,
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
    size_t col_start_x = row * num_individuals;
    size_t col_start_y = col * num_individuals;

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

        float kendall_corr = (concordant - discordant) / sqrt((concordant + discordant + ties_x) *
                                                              (concordant + discordant + ties_y));

        // printf("linear ix: %f, row: %f, col: %f, h: %f, l: %f, corr result: %f \n", lin_ix_f,
        // row,
        //        col, h, l, kendall_corr);

        results[lin_ix] = sin(M_PI / 2 * kendall_corr);
    }
}

// Kendall correlation computation for pairs of markers compressed in .bed format
// without leading magic numbers.
void cu_bed_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                     float *results)
{
    // this is ceil
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t a_bytes = col_len_bytes * num_markers;
    unsigned char *gpu_a;
    float *gpu_results;
    int threads_per_block = NUMTHREADS;
    size_t output_length = num_markers * (num_markers - 1) / 2;
    size_t output_bytes = output_length * sizeof(float);

    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_a, a_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_a, a, a_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_results, output_bytes));

    printf("starting corr kernels \n");

    cu_bed_marker_corr_npn<<<blocks_per_grid, threads_per_block>>>(
        gpu_a, num_markers, num_individuals, col_len_bytes, gpu_results);
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(results, gpu_results, output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_a));
    HANDLE_ERROR(cudaFree(gpu_results));
}

__device__ void unpack_bed_byte(const char b, float *dest)
{
    // TODO: make sure that the bytes are packed from the front,
    // i.e. that the order is most significant -> least significant bits
    for (size_t i = 0; i < 4; i++) {
        dest[i] = bed_lut_a[b + i];
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// compression format is expected to be col-major .bed without NaN.
// that means a should have dimensions ceil(num_individuals / 4) * num_markers
__global__ void cu_bed_marker_corr_npn(const unsigned char *a, const size_t num_markers,
                                       const size_t num_individuals, const size_t col_len_bytes,
                                       float *results)
{
    size_t tix = threadIdx.x;

    // printf("starting corr kernel with id %d \n", tix);

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
    size_t col_start_x = row * col_len_bytes;
    size_t col_start_y = col * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    float bed_vals_x[4] = {0.0};
    float bed_vals_y[4] = {0.0};
    // TODO: it seems stupid to jump in memory, sequential reads are probably more efficient.
    // should have ++ increment and adjust the start.
    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS) {
        // TODO: make sure that unpacking happens in correct order
        printf("block [x: %f; y: %f] thread %d: val of a at %llu: %u \n", col, row, tix,
               col_start_x + i, a[col_start_x + i]);
        printf("block [x: %f; y: %f] thread %d: unpacking byte at x: %llu \n", col, row, tix,
               col_start_x + i);
        unpack_bed_byte(a[col_start_x + i], bed_vals_x);
        printf("block [x: %f; y: %f] thread %d: unpacking byte at y: %llu \n", col, row, tix,
               col_start_y + i);
        unpack_bed_byte(a[col_start_y + i], bed_vals_y);

        for (size_t j = 0; j < 4; j++) {
            if ((i * 4 + j) < num_individuals) {
                size_t comp_ix = (size_t)((3 * bed_vals_x[j] + bed_vals_y[j]));
                thread_sum[comp_ix] += 1.f;
            }
        }
    }

    for (size_t i = 0; i < 9; i++) {
        thread_sums[tix][i] = thread_sum[i];
    }

    // printf("block [x: %f; y: %f] thread %d: finished.", row, col, tix);

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0) {
        // printf("block [x: %f; y: %f]: making single sum", row, col);

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

        float kendall_corr = (concordant - discordant) / sqrt((concordant + discordant + ties_x) *
                                                              (concordant + discordant + ties_y));

        // printf("linear ix: %f, row: %f, col: %f, h: %f, l: %f, corr result: %f \n", lin_ix_f,
        // row,
        //    col, h, l, kendall_corr);

        results[lin_ix] = sin(M_PI / 2 * kendall_corr);
    }
}