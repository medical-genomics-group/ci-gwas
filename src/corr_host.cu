#include <iomanip>
#include <iostream>
#include <math.h>
#include <mps/bed_lut.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/gpuerrors.h>
#include <vector>

void cu_marker_corr_npn_batched(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t row_set_size, // const number of markers
                               // i.e. half the total num markers in a batch
    float *marker_corrs)
{
    size_t num_full_stripes = num_markers / row_set_size;
    size_t last_stripe_width = num_markers % row_set_size;
    bool small_stripe = last_stripe_width;
    size_t num_stripes_total = num_full_stripes + small_stripe;

    size_t num_regular_batches = num_markers / row_set_size - 1;
    size_t ncols_small_batch = num_markers % row_set_size;
    bool small_batch = (ncols_small_batch > 0);

    size_t num_batches = num_regular_batches + small_batch;

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char); // this is ceil
    size_t row_set_bytes = col_len_bytes * row_set_size;

    unsigned char *gpu_marker_vals_row;
    unsigned char *gpu_marker_vals_col;

    float *gpu_marker_corrs;

    size_t batch_output_len = row_set_size * row_set_size;

    // allocate tmp space in host memory
    // batch results are put here before they are put in the full output vec
    // in the correct order
    std::vector<float> marker_corrs_tmp(batch_output_len, 0.0);

    size_t batch_output_bytes = batch_output_len * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // allocate space for marker subsets
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_row, row_set_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_col, row_set_bytes));

    // markers vs markers
    // allocate space for correlation results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, batch_output_bytes));

    // copy row marker data for first stripe to device
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_vals_row,
                   &marker_vals[0],
                   row_set_bytes,
                   cudaMemcpyHostToDevice));

    size_t batch_result_start_host = 0;
    // Iterate through stripes
    // TODO: put loop body in new function, this is hard to read
    for (size_t stripe_ix = 0; stripe_ix < num_stripes_total; ++stripe_ix)
    {
        size_t stripe_width;
        if (small_stripe && stripe_ix == (num_stripes_total - 1))
        {
            stripe_width = last_stripe_width;
        }
        else
        {
            stripe_width = row_set_size;
        }

        size_t stripe_first_row_ix = row_set_size * stripe_ix;
        size_t batch_data_ix = row_set_bytes * (num_regular_batches + 1);         // this is the index of the
                                                                                  // first byte of the first
                                                                                  // column of the first
                                                                                  // (rightmost) batch
                                                                                  // in each stripe
        size_t batch_first_col_ix = row_set_size * (num_regular_batches + 1) - 1; // ix 0 is the first col
                                                                                  // off the diagonal
        size_t batch_data_bytes;
        size_t batch_num_cols;
        size_t batch_num_corrs;
        size_t batch_corrs_bytes;

        if (small_batch)
        {
            batch_num_cols = ncols_small_batch;
        }
        else
        {
            batch_num_cols = row_set_size;
        }

        // Iterate through batches in stripe
        for (size_t batch_ix = 0; batch_ix < num_batches; ++batch_ix)
        {
            batch_data_bytes = batch_num_cols * col_len_bytes;
            batch_num_corrs = batch_num_cols * stripe_width;
            batch_corrs_bytes = batch_num_corrs * sizeof(float);
            dim3 num_blocks(batch_num_cols, stripe_width);

            // copy col marker data to device
            HANDLE_ERROR(
                cudaMemcpy(gpu_marker_vals_col,
                           &marker_vals[batch_data_ix],
                           batch_data_bytes,
                           cudaMemcpyHostToDevice));

            // compute marker correlations
            bed_marker_corr_kendall_npn_batched<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                gpu_marker_vals_col,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs);
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(
                cudaMemcpy(
                    marker_corrs_tmp.data(),
                    gpu_marker_corrs,
                    batch_corrs_bytes,
                    cudaMemcpyDeviceToHost));

            // put correlations from tmp into right place
            size_t rix, cix, lin_ix;
            for (size_t ix = 0; ix < batch_num_corrs; ++ix)
            {
                rix = stripe_first_row_ix + (ix / batch_num_cols);
                cix = batch_first_col_ix + (ix % batch_num_cols);
                lin_ix = (rix * (2 * num_markers - rix - 1) / 2) + (cix - rix); // r (2p - r - 1) / 2 is the first ix of
                                                                                // of the row
                marker_corrs[lin_ix] = marker_corrs_tmp[ix];
            }

            batch_data_ix -= row_set_bytes;
            batch_result_start_host += batch_num_corrs;
            batch_num_cols = row_set_size;
            batch_first_col_ix -= row_set_size;
        }

        // process row markers vs row markers
        if (stripe_width > 1)
        {
            batch_num_corrs = stripe_width * (stripe_width - 1) / 2;
            batch_corrs_bytes = batch_num_corrs * sizeof(float);
            dim3 num_blocks(batch_num_corrs);

            bed_marker_corr_kendall_npn_batched_row<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                stripe_width,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs);
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(
                cudaMemcpy(
                    marker_corrs_tmp.data(),
                    gpu_marker_corrs,
                    batch_corrs_bytes,
                    cudaMemcpyDeviceToHost));

            // put correlations from tmp into right place
            // marker_corrs is upper triangular, batch is upper triangular
            size_t rix, cix, glob_lin_ix;
            size_t loc_lin_ix = 0;
            for (size_t r = 0; r < (stripe_width - 1); ++r)
            {
                rix = stripe_first_row_ix + r;
                for (size_t c = r; c < (stripe_width - 1); ++c)
                {
                    cix = stripe_first_row_ix + c;
                    glob_lin_ix = (rix * (2 * num_markers - rix - 1) / 2) + (cix - rix);
                    marker_corrs[glob_lin_ix] = marker_corrs_tmp[loc_lin_ix];
                    ++loc_lin_ix;
                }
            }

            batch_result_start_host += batch_num_corrs;
        }

        // next stripe is going to have one batch less
        --num_batches;
        // swap gpu_marker_vals_row and gpu_marker_vals_col, bc the last col
        // batch is the next row batch
        unsigned char *tmp = gpu_marker_vals_col;
        gpu_marker_vals_col = gpu_marker_vals_row;
        gpu_marker_vals_row = tmp;
    }

    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
}

void cu_corr_npn_batched(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    const size_t batch_stripe_width, // const number of markers
                                     // i.e. half the total num markers in a batch
    float *marker_corrs,
    float *marker_phen_corrs,
    float *phen_corrs)
{
    size_t num_full_stripes = num_markers / batch_stripe_width;
    size_t last_stripe_width = num_markers % batch_stripe_width;
    bool small_stripe = last_stripe_width;
    size_t num_stripes_total = num_full_stripes + small_stripe;

    size_t num_regular_batches = num_markers / batch_stripe_width - 1; // 1) this checks out because we want
                                                                       //    ((nm - 1) - (bsw - 1)) / bsw =
                                                                       //    (nm - bsw) / bsw = nm / bsw - 1
                                                                       // 2) this is only for the first stripe
                                                                       //    the next stripes will always have
                                                                       //    one less.
    size_t ncols_small_batch = num_markers % batch_stripe_width;
    bool small_batch = (ncols_small_batch > 0);

    size_t num_batches = num_regular_batches + small_batch; // this is only for the first stripe.
                                                            // but it is adjusted later in the loop.

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char); // this is ceil
    size_t batch_marker_vals_bytes = col_len_bytes * batch_stripe_width;
    size_t phen_vals_bytes = num_phen * num_individuals * sizeof(float); // We load all the phenotype data,
                                                                         // no batching here

    unsigned char *gpu_marker_vals_row;
    unsigned char *gpu_marker_vals_col;
    float *gpu_phen_vals;

    float *gpu_marker_corrs;
    float *gpu_marker_phen_corrs;
    float *gpu_phen_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    size_t batch_marker_output_length = batch_stripe_width * batch_stripe_width;
    size_t batch_marker_phen_output_length = batch_stripe_width * num_phen;
    size_t phen_output_length = num_phen * (num_phen - 1) / 2;

    // allocate tmp space in host memory
    // batch results are put here before they are put in the full output vec
    // in the correct order
    std::vector<float> marker_corrs_tmp(batch_marker_output_length, 0.0);
    std::vector<float> marker_phen_corrs_tmp(batch_marker_phen_output_length, 0.0);

    size_t batch_marker_output_bytes = batch_marker_output_length * sizeof(float);
    size_t batch_marker_phen_output_bytes = batch_marker_phen_output_length * sizeof(float);
    size_t phen_output_bytes = phen_output_length * sizeof(float);
    size_t marker_stats_bytes = num_markers * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // allocate space for marker subsets
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_row, batch_marker_vals_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_col, batch_marker_vals_bytes));

    // copy marker stats over, assume that there is always space to fit them in
    // without partitioning
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice));

    // allocate and copy over phenotype data
    HANDLE_ERROR(cudaMalloc(&gpu_phen_vals, phen_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_phen_vals, phen_vals, phen_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_phen_corrs, batch_marker_phen_output_bytes));

    // markers vs markers
    // allocate space for correlation results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, batch_marker_output_bytes));

    // copy row marker data for first stripe to device
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_vals_row,
                   &marker_vals[0],
                   batch_marker_vals_bytes,
                   cudaMemcpyHostToDevice));

    size_t batch_result_start_host = 0;
    // Iterate through stripes
    // TODO: put loop body in new function, this is hard to read
    for (size_t stripe_ix = 0; stripe_ix < num_stripes_total; ++stripe_ix)
    {
        size_t stripe_width;
        if (small_stripe && stripe_ix == (num_stripes_total - 1))
        {
            stripe_width = last_stripe_width;
        }
        else
        {
            stripe_width = batch_stripe_width;
        }

        size_t stripe_first_row_ix = batch_stripe_width * stripe_ix;
        size_t num_marker_phen_corrs = stripe_width * num_phen;
        size_t marker_phen_corrs_bytes = num_marker_phen_corrs * sizeof(float);
        size_t batch_data_ix = batch_marker_vals_bytes * (num_regular_batches + 1);     // this is the index of the
                                                                                        // first byte of the first
                                                                                        // column of the first
                                                                                        // (rightmost) batch
                                                                                        // in each stripe
        size_t batch_first_col_ix = batch_stripe_width * (num_regular_batches + 1) - 1; // ix 0 is the first col
                                                                                        // off the diagonal
        size_t batch_data_bytes;
        size_t batch_num_cols;
        size_t batch_num_corrs;
        size_t batch_corrs_bytes;

        if (small_batch)
        {
            batch_num_cols = ncols_small_batch;
        }
        else
        {
            batch_num_cols = batch_stripe_width;
        }

        // Iterate through batches in stripe
        for (size_t batch_ix = 0; batch_ix < num_batches; ++batch_ix)
        {
            batch_data_bytes = batch_num_cols * col_len_bytes;
            batch_num_corrs = batch_num_cols * stripe_width;
            batch_corrs_bytes = batch_num_corrs * sizeof(float);
            dim3 num_blocks(batch_num_cols, stripe_width);

            // copy col marker data to device
            HANDLE_ERROR(
                cudaMemcpy(gpu_marker_vals_col,
                           &marker_vals[batch_data_ix],
                           batch_data_bytes,
                           cudaMemcpyHostToDevice));

            // compute marker correlations
            bed_marker_corr_kendall_npn_batched<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                gpu_marker_vals_col,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs);
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(
                cudaMemcpy(
                    marker_corrs_tmp.data(),
                    gpu_marker_corrs,
                    batch_corrs_bytes,
                    cudaMemcpyDeviceToHost));

            // put correlations from tmp into right place
            size_t rix, cix, lin_ix;
            for (size_t ix = 0; ix < batch_num_corrs; ++ix)
            {
                rix = stripe_first_row_ix + (ix / batch_num_cols);
                cix = batch_first_col_ix + (ix % batch_num_cols);
                lin_ix = (rix * (2 * num_markers - rix - 1) / 2) + (cix - rix); // r (2p - r - 1) / 2 is the first ix of
                                                                                // of the row
                marker_corrs[lin_ix] = marker_corrs_tmp[ix];
            }

            batch_data_ix -= batch_marker_vals_bytes;
            batch_result_start_host += batch_num_corrs;
            batch_num_cols = batch_stripe_width;
            batch_first_col_ix -= batch_stripe_width;
        }

        // process row markers vs row markers
        if (stripe_width > 1)
        {
            batch_num_corrs = stripe_width * (stripe_width - 1) / 2;
            batch_corrs_bytes = batch_num_corrs * sizeof(float);
            dim3 num_blocks(batch_num_corrs);

            bed_marker_corr_kendall_npn_batched_row<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                stripe_width,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs);
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(
                cudaMemcpy(
                    marker_corrs_tmp.data(),
                    gpu_marker_corrs,
                    batch_corrs_bytes,
                    cudaMemcpyDeviceToHost));

            // put correlations from tmp into right place
            // marker_corrs is upper triangular, batch is upper triangular
            size_t rix, cix, glob_lin_ix;
            size_t loc_lin_ix = 0;
            for (size_t r = 0; r < (stripe_width - 1); ++r)
            {
                rix = stripe_first_row_ix + r;
                for (size_t c = r; c < (stripe_width - 1); ++c)
                {
                    cix = stripe_first_row_ix + c;
                    glob_lin_ix = (rix * (2 * num_markers - rix - 1) / 2) + (cix - rix);
                    marker_corrs[glob_lin_ix] = marker_corrs_tmp[loc_lin_ix];
                    ++loc_lin_ix;
                }
            }

            batch_result_start_host += batch_num_corrs;
        }

        // compute row markers vs phenotype data
        // TODO: this should be Kendall instead of pearson.
        dim3 num_blocks(num_marker_phen_corrs);
        marker_phen_corr_pearson<<<num_blocks, threads_per_block>>>(
            gpu_marker_vals_row,
            gpu_phen_vals,
            stripe_width,
            num_individuals,
            num_phen,
            col_len_bytes,
            &gpu_marker_mean[stripe_first_row_ix],
            &gpu_marker_std[stripe_first_row_ix],
            gpu_marker_phen_corrs);
        CudaCheckError();

        HANDLE_ERROR(
            cudaMemcpy(
                marker_phen_corrs_tmp.data(),
                gpu_marker_phen_corrs,
                marker_phen_corrs_bytes,
                cudaMemcpyDeviceToHost));

        // sort
        // marker_phen_corrs is rectangular (p x num_phen), batch rectangular (stripe_width x num_phen)
        size_t rix, cix;
        for (size_t ix = 0; ix < num_marker_phen_corrs; ++ix)
        {
            rix = stripe_first_row_ix + (ix / num_phen);
            cix = ix % num_phen;
            marker_phen_corrs[cix + (rix * num_phen)] = marker_phen_corrs_tmp[ix];
        }

        // next stripe is going to have one batch less
        --num_batches;
        // swap gpu_marker_vals_row and gpu_marker_vals_col, bc the last col
        // batch is the next row batch
        // std::cerr << "swapping col and row pointers" << std::endl;
        unsigned char *tmp = gpu_marker_vals_col;
        gpu_marker_vals_col = gpu_marker_vals_row;
        gpu_marker_vals_row = tmp;
        // std::cerr << "finished stripe" << std::endl;
    }

    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals_row));
    HANDLE_ERROR(cudaFree(gpu_marker_vals_col));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));

    // phen vs phen
    size_t blocks_per_grid = phen_output_length;
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
