#include <math.h>

#include <mps/gpuerrors.h>
#include <mps/bed_lut.h>
#include <mps/corr_compressed.h>

// Compute correlations between markers, markers and phenotypes, and between phenotypes.
// Markers are expected to be in compressed .bed format, with NaNs removed and no leading magic
// numbers.
void cu_corr_npn(const unsigned char *marker_vals, const float *phen_vals, const size_t num_markers,
                 const size_t num_individuals, const size_t num_phen, const float *marker_mean,
                 const float *marker_std, float *marker_corrs, float *marker_phen_corrs,
                 float *phen_corrs)
{
    // this is ceil
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t marker_vals_bytes = col_len_bytes * num_markers;
    size_t phen_vals_bytes = num_phen * num_individuals * sizeof(float);

    unsigned char *gpu_marker_vals;
    float *gpu_phen_vals;

    float *gpu_marker_corrs;
    float *gpu_marker_phen_corrs;
    float *gpu_phen_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    size_t marker_output_length = num_markers * (num_markers - 1) / 2;
    size_t marker_phen_output_length = num_markers * num_phen;
    size_t phen_output_length = num_phen * (num_phen - 1) / 2;

    size_t marker_output_bytes = marker_output_length * sizeof(float);
    size_t marker_phen_output_bytes = marker_phen_output_length * sizeof(float);
    size_t phen_output_bytes = phen_output_length * sizeof(float);
    size_t marker_stats_bytes = num_markers * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // markers vs markers
    // TODO: see if proper blocks give any performace increase
    int blocks_per_grid = marker_output_length;
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, marker_vals_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, marker_output_bytes));

    bed_marker_corr_kendall_npn<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals, num_markers, num_individuals, col_len_bytes, gpu_marker_corrs);
    CudaCheckError();

    HANDLE_ERROR(
        cudaMemcpy(marker_corrs, gpu_marker_corrs, marker_output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_marker_corrs));

    // markers vs phen
    blocks_per_grid = marker_phen_output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_phen_vals, phen_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_phen_vals, phen_vals, phen_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_phen_corrs, marker_phen_output_bytes));

    marker_phen_corr_pearson<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals, gpu_phen_vals, num_markers, num_individuals, num_phen, col_len_bytes,
        gpu_marker_mean, gpu_marker_std, gpu_marker_phen_corrs);
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(marker_phen_corrs, gpu_marker_phen_corrs, marker_phen_output_bytes,
                            cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_marker_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));

    // phen vs phen
    blocks_per_grid = phen_output_length;
    HANDLE_ERROR(cudaMalloc(&gpu_phen_corrs, phen_output_bytes));

    phen_corr_pearson<<<blocks_per_grid, threads_per_block>>>(gpu_phen_vals, num_individuals,
                                                              num_phen, gpu_phen_corrs);
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(phen_corrs, gpu_phen_corrs, phen_output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_phen_vals));
}

__global__ void marker_phen_corr_pearson(const unsigned char *marker_vals, const float *phen_vals,
                                         const size_t num_markers, const size_t num_individuals,
                                         const size_t num_phen, const size_t col_len_bytes,
                                         const float *marker_mean, const float *marker_std,
                                         float *results)
{
    size_t tix = threadIdx.x;
    size_t lin_ix = blockIdx.x;
    size_t mv_ix = lin_ix / num_phen;
    size_t phen_ix = lin_ix - (num_phen * mv_ix);
    size_t mv_start_ix = mv_ix * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    __shared__ float thread_sums_mv_phen[NUMTHREADS];
    __shared__ float thread_sums_phen[NUMTHREADS];
    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS) {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++) {
            float mv_val = bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];

	    thread_sum_mv_phen += (mv_val * phen_val);
            thread_sum_phen += phen_val;
        }
    }

    thread_sums_mv_phen[tix] = thread_sum_mv_phen;
    thread_sums_phen[tix] = thread_sum_phen;

    __syncthreads();
    if (tix == 0) {
        float s_mv_phen = 0.0;
        float s_phen = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++) {
            s_mv_phen += thread_sums_mv_phen[i];
            s_phen += thread_sums_phen[i];
        }

        results[lin_ix] = (s_mv_phen - marker_mean[mv_ix] * s_phen) /
                          ((float)(num_individuals) * marker_std[mv_ix]);
    }
}

// Compute Pearson's r between pairs of standardized phenotype vectors.
__global__ void phen_corr_pearson(const float *phen_vals, const size_t num_individuals,
                                  const size_t num_phen, float *results)
{
    size_t tix = threadIdx.x;

    // convert linear indices into correlation matrix into (row, col) ix
    size_t lin_ix = blockIdx.x;
    float lin_ix_f = (size_t)lin_ix;
    float l = num_phen - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - lin_ix_f);
    float row = std::floor((-b + sqrt(b * b + 4 * c)) / -2.0) + 1.0;
    float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
    // offset of 1 because we don't compute the diagonal
    float col = lin_ix_f - h + row + 1;
    size_t col_start_a = row * num_individuals;
    size_t col_start_b = col * num_individuals;

    float thread_sum = 0.0;
    __shared__ float thread_sums[NUMTHREADS];
    for (size_t i = tix; i < num_individuals; i += NUMTHREADS) {
        float val_a = phen_vals[(col_start_a + i)];
        float val_b = phen_vals[(col_start_b + i)];
        thread_sum += val_a * val_b;
    }

    thread_sums[tix] = thread_sum;

    __syncthreads();
    if (tix == 0) {
        float s = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++) {
            s += thread_sums[i];
        }
        results[lin_ix] = s / (float)(num_individuals);
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void bed_marker_corr_kendall_npn(const unsigned char *marker_vals,
                                            const size_t num_markers, const size_t num_individuals,
                                            const size_t col_len_bytes, float *results)
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
        size_t aix = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t bix = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++) {
            float val_a = bed_lut_a[(aix + j)];
            float val_b = bed_lut_a[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += 1.f;
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