#include <math.h>

#include <mps/gpuerrors.h>
#include <mps/bed_lut.h>
#include <mps/corr_compressed.h>

void cu_corr_npn_batched(const unsigned char *marker_vals,
                       const float *phen_vals,
                       const size_t num_markers,
                       const size_t num_individuals,
                       const size_t num_phen,
                       const float *marker_mean,
                       const float *marker_std,
                       // the constant number of markers, 
                       // i.e. half the total number of markers in a batch
                       const size_t const_batch_size,
                       float *marker_corrs,
                       float *marker_phen_corrs,
                       float *phen_corrs)
{
    size_t full_width_stripe_length = num_markers - const_batch_size;
    size_t num_regular_batches = full_width_stripe_length / const_batch_size;
    size_t ncols_small_batch = full_width_stripe_length % const_batch_size;
    bool small_batch = (ncols_small_batch == 0);
    // this is ceil
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t marker_vals_bytes = col_len_bytes * const_batch_size;
    size_t phen_vals_bytes = num_phen * num_individuals * sizeof(float);

    unsigned char *gpu_marker_vals;
    float *gpu_phen_vals;

    float *gpu_marker_corrs;
    float *gpu_marker_phen_corrs;
    float *gpu_phen_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    // batch output size
    size_t marker_output_length = const_batch_size * const_batch_size;
    size_t marker_phen_output_length = const_batch_size * num_phen;
    size_t phen_output_length = num_phen * (num_phen - 1) / 2;

    size_t marker_output_bytes = marker_output_length * sizeof(float);
    size_t marker_phen_output_bytes = marker_phen_output_length * sizeof(float);
    size_t phen_output_bytes = phen_output_length * sizeof(float);
    size_t marker_stats_bytes = const_batch_size * sizeof(float);

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

// Compute pearson correlations between markers.
// Markers are expected to be in compressed .bed format, with NaNs removed
// and no leading magic numbers.
void cu_marker_corr_pearson(const unsigned char *marker_vals,
                            const size_t num_markers,
                            const size_t num_individuals,
                            const float *marker_mean,
                            const float *marker_std,
                            float *marker_corrs)
{
    printf("Starting cu_marer_corr_pearson \n");

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t marker_vals_bytes = col_len_bytes * num_markers;

    unsigned char *gpu_marker_vals;
    float *gpu_marker_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    size_t marker_output_length = num_markers * (num_markers - 1) / 2;
    size_t marker_output_bytes = marker_output_length * sizeof(float);
    size_t marker_stats_bytes = num_markers * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // markers vs markers
    int blocks_per_grid = marker_output_length;
    
    
    // allocate device memory and copy over data
    // marker values
    printf("Allocating space for marker data \n");
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, marker_vals_bytes)); 
    printf("Copying marker data to device \n");
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice));

    // space for results
    printf("Allocating space for results \n");
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, marker_output_bytes));

    // marker means
    printf("Allocating space for marker means \n");
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    printf("Copying marker means to device \n");
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice));

    // marker stds
    printf("Allocating space for marker stds \n");
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes)); 
    printf("Copying marker stds to device \n");
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice));

    printf("Calling kernels \n");

    // compute correlations
    marker_corr_pearson<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals, num_markers, num_individuals,
        col_len_bytes, gpu_marker_mean, gpu_marker_std, gpu_marker_corrs);
    CudaCheckError();
    
    printf("Copying results to host \n");
    // copy results to host
    HANDLE_ERROR(cudaMemcpy(marker_corrs, gpu_marker_corrs, marker_output_bytes, cudaMemcpyDeviceToHost));
    
    
    printf("Freeing device memory \n");
    // free allocated device memory
    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));
}

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

// Compute row and column indices from the linear index into
// upper triangular matrix.
__device__ void row_col_ix_from_linear_ix(const size_t lin_ix,
                                          const size_t num_rows,
                                          size_t *row_ix,
                                          size_t *col_ix)
{
    float l = num_rows - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - (float)lin_ix);
    float row = std::floor((-b + sqrt(b * b + 4 * c)) / - 2.0) + 1.0;
    float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
    // offset of 1 because we do not include the diagonal
    float col = (float)lin_ix - h + row + 1;
    *row_ix = (size_t)row;
    *col_ix = (size_t)col;
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
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
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

// Compute Pearson's r between a pair of standardized phenotype vectors.
__global__ void phen_corr_pearson(const float *phen_vals, const size_t num_individuals,
                                  const size_t num_phen, float *results)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_phen, &row, &col);
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
        results[blockIdx.x] = s / (float)(num_individuals);
    }
}

// Compute Pearson's r between a pair of marker vectors.
__global__ void marker_corr_pearson(const unsigned char *marker_vals,
                                    const size_t num_markers,
                                    const size_t num_individuals,
                                    const size_t col_len_bytes,
                                    const float *marker_mean,
                                    const float *marker_std,
                                    float *results)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum_mvp = 0.0;
    __shared__ float thread_sums_mvp[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS) {
        size_t mv_byte_ix_a = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t mv_byte_ix_b = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++) {
                float mv_a = gpu_bed_lut_a[(mv_byte_ix_a + j)];
                float mv_b = gpu_bed_lut_a[(mv_byte_ix_b + j)];
                thread_sum_mvp += mv_a * mv_b;
        }
    }

    thread_sums_mvp[tix] = thread_sum_mvp;

    __syncthreads();
    //printf("[tix: %lu]: syncing... \n", tix);
    if (tix == 0) {
        float s_mvps = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++) {
            s_mvps += thread_sums_mvp[i];
        }

        float num = (s_mvps / (float)num_individuals) - (marker_mean[row] * marker_mean[col]);
        float denom = marker_std[row] * marker_std[col];

        results[blockIdx.x] = num / denom;
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
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
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
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
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

        results[blockIdx.x] = sin(M_PI / 2 * kendall_corr);
    }
}
