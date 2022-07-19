#include <math.h>
#include <mps/bed_lut.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/gpuerrors.h>

void cu_corr_npn_batched(const unsigned char *marker_vals,
                         const float *phen_vals,
                         const size_t num_markers,
                         const size_t num_individuals,
                         const size_t num_phen,
                         const float *marker_mean,
                         const float *marker_std,
                         // the constant number of markers,
                         // i.e. half the total number of markers in a batch
                         const size_t row_width,
                         float *marker_corrs,
                         float *marker_phen_corrs,
                         float *phen_corrs)
{
    size_t num_full_stripes = (num_markers - 1) / row_width;
    size_t nrows_small_stripe = (num_markers - 1) % row_width;
    bool small_stripe = nrows_small_stripe;
    size_t num_stripes = num_full_stripes + small_stripe;

    size_t num_regular_batches = num_full_stripes / row_width - 1;
    size_t ncols_small_batch = num_full_stripes % row_width;
    bool small_batch = (ncols_small_batch == 0);
    size_t num_batches = num_regular_batches + small_batch;

    // this is ceil
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t marker_vals_bytes = col_len_bytes * row_width;
    size_t phen_vals_bytes = num_phen * num_individuals * sizeof(float);

    unsigned char *gpu_marker_vals_row;
    unsigned char *gpu_marker_vals_col;
    float *gpu_phen_vals;

    float *gpu_marker_corrs;
    float *gpu_marker_phen_corrs;
    float *gpu_phen_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    // batch output size
    size_t marker_output_length = row_width * row_width;
    size_t marker_phen_output_length = row_width * num_phen;
    size_t phen_output_length = num_phen * (num_phen - 1) / 2;

    size_t marker_output_bytes = marker_output_length * sizeof(float);
    size_t marker_phen_output_bytes = marker_phen_output_length * sizeof(float);
    size_t phen_output_bytes = phen_output_length * sizeof(float);
    // for now we don't partition the stats
    size_t marker_stats_bytes = num_markers * sizeof(float);

    int threads_per_block = NUMTHREADS;

    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_row, marker_vals_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_col, marker_vals_bytes));

    // copy marker stats over, assume that there is always space to fit them in
    // without partitioning
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice));

    // markers vs markers
    for (size_t stripe_ix = 0; stripe_ix < num_stripes; ++stripe_ix)
    {
        size_t stripe_marker_bytes, stripe_width;
        if (small_stripe && stripe_ix == (num_stripes - 1))
        {
            stripe_marker_bytes = (nrows_small_stripe * col_len_bytes);
            stripe_width = nrows_small_stripe;
        }
        else
        {
            stripe_marker_bytes = marker_vals_bytes;
            stripe_width = row_width;
        }

        // copy row marker data to device
        HANDLE_ERROR(
            cudaMemcpy(gpu_marker_vals_row,
                       marker_vals[stripe_ix * col_len_bytes],
                       stripe_marker_bytes,
                       cudaMemcpyHostToDevice));

        size_t batch_data_ix = 0;
        size_t batchbytes;
        for (size_t batch_ix = 0; batch_ix < num_batches; ++batch_ix)
        {
            size_t nbytes;
            dim3 num_blocks;

            if (small_batch && batch_ix == 0)
            {
                batchbytes = ncols_small_batch * col_len_bytes;
                ;
                num_blocks = dim3(ncols_small_batch, stripe_width);
            }
            else
            {
                batchbytes = row_width * col_len_bytes;
                num_blocks = dim3(row_width, stripe_width);
            }

            // copy col marker data to device
            HANDLE_ERROR(
                cudaMemcpy(gpu_marker_vals_col,
                           marker_vals[batch_data_ix],
                           batchbytes,
                           cudaMemcpyHostToDevice));

            // at this point I have marker data for both row and cols
            // and all stats, so I can call a kernel here.
            // I am computing a rectangular matrix thing, so can use 2D blocks

            batch_data_ix += batchbytes;
        }
        // don't forget to process that row elements against themselves!

        // next stripe is going to have one batch less
        --num_batches;
        // swap gpu_marker_vals_row and gpu_marker_vals_col, bc the last col
        // batch is the next row batch
    }

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