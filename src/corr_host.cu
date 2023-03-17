#include <math.h>
#include <mps/bed_lut_gpu.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/gpuerrors.h>
#include <mps/io.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

using byte = unsigned char;

std::vector<float> cal_mcorrk(
    const std::vector<unsigned char> &bed_vals,
    const BedDims dim,
    const size_t num_markers,
    const size_t device_mem_gb
)
{
    size_t num_individuals = dim.get_num_samples();
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    // mem required for non-batched processing
    size_t nbytes_marker_data = num_individuals / 4 * num_markers;
    size_t nbytes_marker_corrs = num_markers * (num_markers - 1) * 2;
    size_t req_mem_bytes = nbytes_marker_data + nbytes_marker_corrs;

    if (req_mem_bytes > device_mem_bytes)
    {
        // figure out batch size
        float b = num_individuals / 4 - 2;
        size_t max_batch_size = (size_t)((-b + std::sqrt(b * b + 8 * device_mem_bytes)) / 4);
        size_t batch_nrows = max_batch_size / 2;

        printf("req mem: %zu, device mem: %zu \n", req_mem_bytes, device_mem_bytes);
        size_t mb_bytes =
            num_individuals / 4 * max_batch_size + max_batch_size * (max_batch_size - 1) * 2;
        printf("max_batch_size: %zu (%zu bytes \n)\n", max_batch_size, mb_bytes);
        printf("Device mem < required mem; Running tiled routine. \n");
        fflush(stdout);

        cu_marker_corr_pearson_npn_batched(
            bed_vals.data(), num_markers, num_individuals, batch_nrows, marker_corr.data()
        );
    }
    else
    {
        printf("Device mem >= required mem; Running non-tiled routine. \n");
        fflush(stdout);
        // compute correlations
        cu_marker_corr_pearson_npn(
            bed_vals.data(), num_markers, num_individuals, marker_corr.data()
        );
    }

    return marker_corr;
}

std::vector<float> cal_mcorrk_banded(
    const std::vector<unsigned char> &bed_vals,
    const BedDims dim,
    const size_t num_markers,
    const size_t corr_width,
    const size_t device_mem_gb
)
{
    size_t num_individuals = dim.get_num_samples();
    float device_mem_bytes = device_mem_gb * std::pow(10, 9);
    size_t upper_bound_batch_size = std::floor((device_mem_bytes) / ((num_individuals + 3) / 4));

    if (upper_bound_batch_size > num_markers)
    {
        upper_bound_batch_size = num_markers;
    }
    if (upper_bound_batch_size < corr_width)
    {
        printf(
            "Maximal batch size (%u) < corr width (%u). Decrease distance threshold or increase "
            "device memory. \n",
            upper_bound_batch_size,
            corr_width
        );
        exit(1);
    }

    // allocate correlation result arrays
    size_t corr_mat_size = corr_width * num_markers;
    std::vector<float> corrs(corr_mat_size, 0.0);

    // compute correlations
    cu_marker_corr_pearson_npn_batched_sparse(
        bed_vals.data(),
        num_markers,
        num_individuals,
        corr_width,
        upper_bound_batch_size,
        corrs.data()
    );

    return corrs;
}

std::vector<float> marker_corr_banded_mat_antidiag_sums(
    const std::vector<float> &marker_corrs, size_t corr_width
)
{
    size_t num_markers = marker_corrs.size() / corr_width;
    size_t num_anti_diagonals = num_markers - 1;
    std::vector<float> sums(num_anti_diagonals, 0.0);

    for (size_t row = 0; row < num_markers - 1; row++)
    {
        for (size_t col = 0; col < corr_width; col++)
        {
            size_t diag_ix = row + col;
            if (diag_ix < num_anti_diagonals)
            {
                sums[row + col] += marker_corrs[row * corr_width + col];
            }
        }
    }

    return sums;
}

std::vector<float> marker_corr_mat_antidiag_sums(const std::vector<float> &marker_corrs)
{
    size_t num_markers = (marker_corrs.size() + 3) / 2;
    size_t row_start_ix = 0;
    size_t row_size = num_markers - 1;
    size_t num_anti_diagonals = 2 * num_markers - 3;
    std::vector<float> sums(num_anti_diagonals, 0.0);

    for (size_t row = 0; row < (num_markers - 1); row++)
    {
        for (size_t update_ix = 0; update_ix < row_size; update_ix++)
        {
            sums[2 * row + update_ix] += marker_corrs[row_start_ix + update_ix];
        }

        row_start_ix += row_size;
        --row_size;
    }

    return sums;
}

void marker_corr_mat_antidiag_sums(const size_t num_markers, const float *marker_corrs, float *sums)
{
    size_t row_start_ix = 0;
    size_t row_size = num_markers - 1;

    for (size_t row = 0; row < (num_markers - 1); row++)
    {
        for (size_t update_ix = 0; update_ix < row_size; update_ix++)
        {
            sums[2 * row + update_ix] += marker_corrs[row_start_ix + update_ix];
        }

        row_start_ix += row_size;
        --row_size;
    }
}

void cu_ix_from_linear(const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix)
{
    size_t *gpu_rix, *gpu_cix;
    HANDLE_ERROR(cudaMalloc(&gpu_rix, sizeof(size_t)));
    HANDLE_ERROR(cudaMalloc(&gpu_cix, sizeof(size_t)));
    ix_from_linear<<<1, 1>>>(lin_ix, num_rows, gpu_rix, gpu_cix);
    CudaCheckError();
    HANDLE_ERROR(cudaMemcpy(row_ix, gpu_rix, sizeof(size_t), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(col_ix, gpu_cix, sizeof(size_t), cudaMemcpyDeviceToHost));
}

void cu_marker_corr_pearson_batched(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const float *marker_mean,
    const float *marker_std,
    const size_t batch_stripe_width,
    float *marker_corrs
)
{
    size_t num_full_stripes = num_markers / batch_stripe_width;
    size_t last_stripe_width = num_markers % batch_stripe_width;
    bool small_stripe = last_stripe_width;
    size_t num_stripes_total = num_full_stripes + small_stripe;

    size_t num_regular_batches = num_markers / batch_stripe_width - 1;
    size_t ncols_small_batch = num_markers % batch_stripe_width;
    bool small_batch = (ncols_small_batch > 0);

    size_t num_batches =
        num_regular_batches + small_batch;  // this is only for the first stripe.
                                            // but it is adjusted later in the loop.

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);  // this is ceil
    size_t batch_marker_vals_bytes = col_len_bytes * batch_stripe_width;

    unsigned char *gpu_marker_vals_row;
    unsigned char *gpu_marker_vals_col;

    float *gpu_marker_corrs;
    float *gpu_marker_mean;
    float *gpu_marker_std;

    size_t batch_marker_output_length = batch_stripe_width * batch_stripe_width;

    // allocate tmp space in host memory
    // batch results are put here before they are put in the full output vec
    // in the correct order
    std::vector<float> marker_corrs_tmp(batch_marker_output_length, 0.0);

    size_t batch_marker_output_bytes = batch_marker_output_length * sizeof(float);
    size_t marker_stats_bytes = num_markers * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // allocate space for marker subsets
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_row, batch_marker_vals_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals_col, batch_marker_vals_bytes));

    // copy marker stats over, assume that there is always space to fit them in
    // without partitioning
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice)
    );

    // markers vs markers
    // allocate space for correlation results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, batch_marker_output_bytes));

    // copy row marker data for first stripe to device
    HANDLE_ERROR(cudaMemcpy(
        gpu_marker_vals_row, &marker_vals[0], batch_marker_vals_bytes, cudaMemcpyHostToDevice
    ));

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
        size_t batch_data_ix = batch_marker_vals_bytes * (num_regular_batches + 1);
        size_t batch_first_col_ix = batch_stripe_width * (num_regular_batches + 1) - 1;
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
            HANDLE_ERROR(cudaMemcpy(
                gpu_marker_vals_col,
                &marker_vals[batch_data_ix],
                batch_data_bytes,
                cudaMemcpyHostToDevice
            ));

            // compute marker correlations
            bed_marker_corr_pearson_batched<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                gpu_marker_vals_col,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                &gpu_marker_mean[stripe_first_row_ix],
                &gpu_marker_std[stripe_first_row_ix],
                &gpu_marker_mean[batch_first_col_ix + 1],
                &gpu_marker_std[batch_first_col_ix + 1],
                gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

            // put correlations from tmp into right place
            size_t rix, cix, lin_ix;
            for (size_t ix = 0; ix < batch_num_corrs; ++ix)
            {
                rix = stripe_first_row_ix + (ix / batch_num_cols);
                cix = batch_first_col_ix + (ix % batch_num_cols);
                lin_ix = (rix * (2 * num_markers - rix - 1) / 2) +
                         (cix - rix);  // r (2p - r - 1) / 2 is the first ix of
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

            bed_marker_corr_pearson_batched_row<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                &gpu_marker_mean[stripe_first_row_ix],
                &gpu_marker_std[stripe_first_row_ix],
                gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

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
        // std::cerr << "swapping col and row pointers" << std::endl;
        unsigned char *tmp = gpu_marker_vals_col;
        gpu_marker_vals_col = gpu_marker_vals_row;
        gpu_marker_vals_row = tmp;
        // std::cerr << "finished stripe" << std::endl;
    }

    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals_row));
    HANDLE_ERROR(cudaFree(gpu_marker_vals_col));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));
}

void cu_marker_corr_pearson_npn_batched(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t row_set_size,  // const number of markers
                                // i.e. half the total num markers in a batch
    float *marker_corrs
)
{
    size_t num_full_stripes = num_markers / row_set_size;
    size_t last_stripe_width = num_markers % row_set_size;
    bool small_stripe = last_stripe_width;
    size_t num_stripes_total = num_full_stripes + small_stripe;

    size_t num_regular_batches = num_markers / row_set_size - 1;
    size_t ncols_small_batch = num_markers % row_set_size;
    bool small_batch = (ncols_small_batch > 0);

    size_t num_batches = num_regular_batches + small_batch;

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);  // this is ceil
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
        cudaMemcpy(gpu_marker_vals_row, &marker_vals[0], row_set_bytes, cudaMemcpyHostToDevice)
    );

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
        size_t batch_data_ix = row_set_bytes * (num_regular_batches + 1);  // this is the index of
                                                                           // the first byte of the
                                                                           // first column of the
                                                                           // first (rightmost)
                                                                           // batch in each stripe
        size_t batch_first_col_ix =
            row_set_size * (num_regular_batches + 1) - 1;  // ix 0 is the first col
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
            HANDLE_ERROR(cudaMemcpy(
                gpu_marker_vals_col,
                &marker_vals[batch_data_ix],
                batch_data_bytes,
                cudaMemcpyHostToDevice
            ));

            // compute marker correlations
            bed_marker_corr_pearson_npn_batched<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                gpu_marker_vals_col,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

            // put correlations from tmp into right place
            size_t rix, cix, lin_ix;
            for (size_t ix = 0; ix < batch_num_corrs; ++ix)
            {
                rix = stripe_first_row_ix + (ix / batch_num_cols);
                cix = batch_first_col_ix + (ix % batch_num_cols);
                lin_ix = (rix * (2 * num_markers - rix - 1) / 2) +
                         (cix - rix);  // r (2p - r - 1) / 2 is the first ix of
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

            bed_marker_corr_pearson_npn_batched_row<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row, stripe_width, num_individuals, col_len_bytes, gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

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
    HANDLE_ERROR(cudaFree(gpu_marker_vals_row));
    HANDLE_ERROR(cudaFree(gpu_marker_vals_col));
}

void cu_corr_pearson_npn_batched(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    const size_t batch_stripe_width,  // const number of markers
                                      // i.e. half the total num markers in a batch
    float *marker_corrs,
    float *marker_phen_corrs,
    float *phen_corrs
)
{
    size_t num_full_stripes = num_markers / batch_stripe_width;
    size_t last_stripe_width = num_markers % batch_stripe_width;
    bool small_stripe = last_stripe_width;
    size_t num_stripes_total = num_full_stripes + small_stripe;

    size_t num_regular_batches =
        num_markers / batch_stripe_width - 1;  // 1) this checks out because we want
                                               //    ((nm - 1) - (bsw - 1)) / bsw =
                                               //    (nm - bsw) / bsw = nm / bsw - 1
                                               // 2) this is only for the first stripe
                                               //    the next stripes will always have
                                               //    one less.
    size_t ncols_small_batch = num_markers % batch_stripe_width;
    bool small_batch = (ncols_small_batch > 0);

    size_t num_batches =
        num_regular_batches + small_batch;  // this is only for the first stripe.
                                            // but it is adjusted later in the loop.

    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);  // this is ceil
    size_t batch_marker_vals_bytes = col_len_bytes * batch_stripe_width;
    size_t phen_vals_bytes =
        num_phen * num_individuals * sizeof(float);  // We load all the phenotype data,
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
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice)
    );

    // allocate and copy over phenotype data
    HANDLE_ERROR(cudaMalloc(&gpu_phen_vals, phen_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_phen_vals, phen_vals, phen_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_phen_corrs, batch_marker_phen_output_bytes));

    // markers vs markers
    // allocate space for correlation results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, batch_marker_output_bytes));

    // copy row marker data for first stripe to device
    HANDLE_ERROR(cudaMemcpy(
        gpu_marker_vals_row, &marker_vals[0], batch_marker_vals_bytes, cudaMemcpyHostToDevice
    ));

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
        size_t batch_data_ix =
            batch_marker_vals_bytes * (num_regular_batches + 1);  // this is the index of the
                                                                  // first byte of the first
                                                                  // column of the first
                                                                  // (rightmost) batch
                                                                  // in each stripe
        size_t batch_first_col_ix =
            batch_stripe_width * (num_regular_batches + 1) - 1;  // ix 0 is the first col
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
            HANDLE_ERROR(cudaMemcpy(
                gpu_marker_vals_col,
                &marker_vals[batch_data_ix],
                batch_data_bytes,
                cudaMemcpyHostToDevice
            ));

            // compute marker correlations
            bed_marker_corr_pearson_npn_batched<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row,
                gpu_marker_vals_col,
                batch_num_cols,
                num_individuals,
                col_len_bytes,
                gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

            // put correlations from tmp into right place
            size_t rix, cix, lin_ix;
            for (size_t ix = 0; ix < batch_num_corrs; ++ix)
            {
                rix = stripe_first_row_ix + (ix / batch_num_cols);
                cix = batch_first_col_ix + (ix % batch_num_cols);
                lin_ix = (rix * (2 * num_markers - rix - 1) / 2) +
                         (cix - rix);  // r (2p - r - 1) / 2 is the first ix of
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

            bed_marker_corr_pearson_npn_batched_row<<<num_blocks, threads_per_block>>>(
                gpu_marker_vals_row, stripe_width, num_individuals, col_len_bytes, gpu_marker_corrs
            );
            CudaCheckError();

            // copy corr results to host
            HANDLE_ERROR(cudaMemcpy(
                marker_corrs_tmp.data(), gpu_marker_corrs, batch_corrs_bytes, cudaMemcpyDeviceToHost
            ));

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
        bed_marker_phen_corr_pearson<<<num_blocks, threads_per_block>>>(
            gpu_marker_vals_row,
            gpu_phen_vals,
            stripe_width,
            num_individuals,
            num_phen,
            col_len_bytes,
            &gpu_marker_mean[stripe_first_row_ix],
            &gpu_marker_std[stripe_first_row_ix],
            gpu_marker_phen_corrs
        );
        CudaCheckError();

        HANDLE_ERROR(cudaMemcpy(
            marker_phen_corrs_tmp.data(),
            gpu_marker_phen_corrs,
            marker_phen_corrs_bytes,
            cudaMemcpyDeviceToHost
        ));

        // sort
        // marker_phen_corrs is rectangular (p x num_phen), batch rectangular (stripe_width x
        // num_phen)
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

    phen_corr_pearson<<<blocks_per_grid, threads_per_block>>>(
        gpu_phen_vals, num_individuals, num_phen, gpu_phen_corrs
    );
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(phen_corrs, gpu_phen_corrs, phen_output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_phen_vals));
}

// Compute pearson correlations between markers.
// Markers are expected to be in compressed .bed format, with NaNs removed
// and no leading magic numbers.
void cu_marker_corr_pearson(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const float *marker_mean,
    const float *marker_std,
    float *marker_corrs
)
{
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
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, marker_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice)
    );

    // space for results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, marker_output_bytes));

    // marker means
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice)
    );

    // marker stds
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice)
    );

    // compute correlations
    bed_marker_corr_pearson<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals,
        num_markers,
        num_individuals,
        col_len_bytes,
        gpu_marker_mean,
        gpu_marker_std,
        gpu_marker_corrs
    );
    CudaCheckError();

    // copy results to host
    HANDLE_ERROR(
        cudaMemcpy(marker_corrs, gpu_marker_corrs, marker_output_bytes, cudaMemcpyDeviceToHost)
    );

    // free allocated device memory
    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));
}

// Compute pearson correlations between markers as sin(pi / 2 * tau-b).
// Markers are expected to be in compressed .bed format, with NaNs removed
// and no leading magic numbers.
void cu_marker_corr_pearson_npn(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    float *marker_corrs
)
{
    size_t col_len_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t marker_vals_bytes = col_len_bytes * num_markers;

    unsigned char *gpu_marker_vals;
    float *gpu_marker_corrs;

    size_t marker_output_length = num_markers * (num_markers - 1) / 2;
    size_t marker_output_bytes = marker_output_length * sizeof(float);

    int threads_per_block = NUMTHREADS;

    // markers vs markers
    int blocks_per_grid = marker_output_length;

    // allocate device memory and copy over data
    // marker values
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, marker_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice)
    );

    // space for results
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, marker_output_bytes));

    // compute correlations
    bed_marker_corr_pearson_npn<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals, num_markers, num_individuals, col_len_bytes, gpu_marker_corrs
    );
    CudaCheckError();

    // copy results to host
    HANDLE_ERROR(
        cudaMemcpy(marker_corrs, gpu_marker_corrs, marker_output_bytes, cudaMemcpyDeviceToHost)
    );

    // free allocated device memory
    HANDLE_ERROR(cudaFree(gpu_marker_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
}

// Compute correlations between markers, markers and phenotypes, and between phenotypes.
void cu_corr_pearson_npn(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    float *marker_corrs,
    float *marker_phen_corrs,
    float *phen_corrs
)
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
    HANDLE_ERROR(cudaMemcpy(gpu_marker_vals, marker_vals, marker_vals_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMalloc(&gpu_marker_corrs, marker_output_bytes));

    bed_marker_corr_pearson_npn_scan<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals, num_markers, num_individuals, col_len_bytes, gpu_marker_corrs
    );
    CudaCheckError();

    HANDLE_ERROR(
        cudaMemcpy(marker_corrs, gpu_marker_corrs, marker_output_bytes, cudaMemcpyDeviceToHost)
    );
    HANDLE_ERROR(cudaFree(gpu_marker_corrs));

    // markers vs phen
    blocks_per_grid = marker_phen_output_length;

    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMalloc(&gpu_phen_vals, phen_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_phen_vals, phen_vals, phen_vals_bytes, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_phen_corrs, marker_phen_output_bytes));

    bed_marker_phen_corr_pearson_scan<<<blocks_per_grid, threads_per_block>>>(
        gpu_marker_vals,
        gpu_phen_vals,
        num_markers,
        num_individuals,
        num_phen,
        col_len_bytes,
        gpu_marker_mean,
        gpu_marker_std,
        gpu_marker_phen_corrs
    );
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(
        marker_phen_corrs, gpu_marker_phen_corrs, marker_phen_output_bytes, cudaMemcpyDeviceToHost
    ));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_marker_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));

    // phen vs phen
    blocks_per_grid = phen_output_length;
    HANDLE_ERROR(cudaMalloc(&gpu_phen_corrs, phen_output_bytes));

    phen_corr_pearson_scan<<<blocks_per_grid, threads_per_block>>>(
        gpu_phen_vals, num_individuals, num_phen, gpu_phen_corrs
    );
    CudaCheckError();

    HANDLE_ERROR(cudaMemcpy(phen_corrs, gpu_phen_corrs, phen_output_bytes, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(gpu_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_phen_vals));
}

void cu_marker_corr_pearson_npn_batched_sparse(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t corr_width,
    const size_t batch_size,
    float *corrs
)
{
    dim3 BLOCKS_PER_GRID;
    dim3 THREADS_PER_BLOCK;

    float *gpu_corrs;
    HANDLE_ERROR(cudaMalloc(&gpu_corrs, corr_width * batch_size * sizeof(float)));

    // copy first marker batch to device
    byte *gpu_marker_vals;
    size_t genotype_col_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t batch_marker_vals_bytes = genotype_col_bytes * batch_size;
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, batch_marker_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(
        gpu_marker_vals, &marker_vals[0], batch_marker_vals_bytes, cudaMemcpyHostToDevice
    ));

    size_t row_in = 0;
    size_t corr_row_len = corr_width;
    size_t last_full_row = num_markers - corr_width - 1;
    size_t max_row_in = batch_size - corr_width;
    size_t batch_start = 0;
    size_t curr_width = corr_width;
    for (size_t row_ix = 0; row_ix < num_markers; ++row_ix)
    {
        if (row_ix % 10000 == 0)
        {
            printf("Processing marker #%u \n", row_ix);
            fflush(stdout);
        }
        size_t row_out = row_ix % batch_size;
        if (curr_width < corr_width)
        {
            // overwrite old results
            HANDLE_ERROR(cudaMemset(
                &gpu_corrs[row_out * corr_row_len + curr_width],
                0,
                (corr_width - curr_width) * sizeof(float)
            ));
        }

        // marker-marker corr
        BLOCKS_PER_GRID = dim3(curr_width, 1, 1);
        THREADS_PER_BLOCK = dim3(NUMTHREADS, 1, 1);
        if (curr_width > 0)
        {
            bed_marker_corr_pearson_npn_sparse_scan<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                gpu_marker_vals,
                num_individuals,
                genotype_col_bytes,
                row_in,
                &gpu_corrs[row_out * corr_row_len]
            );
            CudaCheckError();
        }

        if (row_out == (batch_size - 1))
        {
            printf("Copying batch results to host\n");
            fflush(stdout);
            // copy batch results to host
            HANDLE_ERROR(cudaMemcpy(
                &corrs[batch_start * (corr_width)],
                gpu_corrs,
                (corr_width)*batch_size * sizeof(float),
                cudaMemcpyDeviceToHost
            ));
            batch_start = row_ix + 1;
            fflush(stdout);
        }

        if (row_ix >= last_full_row)
        {
            curr_width -= 1;
            max_row_in += 1;
        }

        if (row_in == (max_row_in - 1))
        {
            printf("At marker #%u;\n", row_ix);
            printf("Loading new marker batch \n");
            fflush(stdout);
            // get new marker data batch
            HANDLE_ERROR(cudaMemcpy(
                gpu_marker_vals,
                &marker_vals[genotype_col_bytes * (row_ix + 1)],
                batch_marker_vals_bytes,
                cudaMemcpyHostToDevice
            ));
            row_in = 0;
        }
        else
        {
            row_in += 1;
        }
    }

    // copy remaining results to host
    if ((num_markers % batch_size) != 0)
    {
        size_t last_batch_size = num_markers - batch_start;
        HANDLE_ERROR(cudaMemcpy(
            &corrs[batch_start * corr_width],
            gpu_corrs,
            corr_width * last_batch_size * sizeof(float),
            cudaMemcpyDeviceToHost
        ));
    }

    HANDLE_ERROR(cudaFree(gpu_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
}

void cu_corr_pearson_npn_batched_sparse(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    const size_t corr_width,
    const size_t batch_size,
    float *corrs
)
{
    /*
    m: num_markers
    m': batch_size
    n: num_individuals
    p: num_phen
    w: corr_width

    memory requirement is 4(np + m'(w + p)) + nm' / 4 bytes
    */

    dim3 BLOCKS_PER_GRID;
    dim3 THREADS_PER_BLOCK;

    assert((batch_size > corr_width) && "error: batch_size <= corr_width");

    // copy all marker stats over
    float *gpu_marker_mean;
    float *gpu_marker_std;
    size_t marker_stats_bytes = num_markers * sizeof(float);
    HANDLE_ERROR(cudaMalloc(&gpu_marker_mean, marker_stats_bytes));
    HANDLE_ERROR(cudaMalloc(&gpu_marker_std, marker_stats_bytes));
    HANDLE_ERROR(
        cudaMemcpy(gpu_marker_mean, marker_mean, marker_stats_bytes, cudaMemcpyHostToDevice)
    );
    HANDLE_ERROR(cudaMemcpy(gpu_marker_std, marker_std, marker_stats_bytes, cudaMemcpyHostToDevice)
    );

    // allocate and copy over phenotype data
    float *gpu_phen_vals;
    size_t phen_vals_bytes = num_phen * num_individuals * sizeof(float);
    HANDLE_ERROR(cudaMalloc(&gpu_phen_vals, phen_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(gpu_phen_vals, phen_vals, phen_vals_bytes, cudaMemcpyHostToDevice));

    // TODO: This could be more than one row, then I could make larger transfers
    // allocate space for correlation results
    float *gpu_corrs;
    HANDLE_ERROR(cudaMalloc(&gpu_corrs, (corr_width + num_phen) * batch_size * sizeof(float)));

    // copy first marker batch to device
    byte *gpu_marker_vals;
    size_t genotype_col_bytes = (num_individuals + 3) / 4 * sizeof(unsigned char);
    size_t batch_marker_vals_bytes = genotype_col_bytes * batch_size;
    HANDLE_ERROR(cudaMalloc(&gpu_marker_vals, batch_marker_vals_bytes));
    HANDLE_ERROR(cudaMemcpy(
        gpu_marker_vals, &marker_vals[0], batch_marker_vals_bytes, cudaMemcpyHostToDevice
    ));

    size_t row_in = 0;
    size_t corr_row_len = corr_width + num_phen;
    size_t last_full_row = num_markers - corr_width - 1;
    size_t max_row_in = batch_size - corr_width;
    size_t batch_start = 0;
    size_t curr_width = corr_width;
    for (size_t row_ix = 0; row_ix < num_markers; ++row_ix)
    {
        if (row_ix % 10000 == 0)
        {
            printf("Processing marker #%u \n", row_ix);
            fflush(stdout);
        }
        size_t row_out = row_ix % batch_size;
        if (curr_width < corr_width)
        {
            // overwrite old results
            HANDLE_ERROR(cudaMemset(
                &gpu_corrs[row_out * corr_row_len + curr_width],
                0,
                (corr_width - curr_width) * sizeof(float)
            ));
        }

        // marker-marker corr
        BLOCKS_PER_GRID = dim3(curr_width, 1, 1);
        THREADS_PER_BLOCK = dim3(NUMTHREADS, 1, 1);
        if (curr_width > 0)
        {
            bed_marker_corr_pearson_npn_sparse_scan<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                gpu_marker_vals,
                num_individuals,
                genotype_col_bytes,
                row_in,
                &gpu_corrs[row_out * corr_row_len]
            );
            CudaCheckError();
        }
        // marker-phen corr
        BLOCKS_PER_GRID = dim3(num_phen, 1, 1);
        THREADS_PER_BLOCK = dim3(NUMTHREADS, 1, 1);
        bed_marker_phen_corr_pearson_sparse_scan<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
            gpu_marker_vals,
            gpu_phen_vals,
            num_individuals,
            num_phen,
            genotype_col_bytes,
            &gpu_marker_mean[row_ix],
            &gpu_marker_std[row_ix],
            corr_width,
            row_in,
            &gpu_corrs[row_out * corr_row_len]
        );
        CudaCheckError();

        if (row_out == (batch_size - 1))
        {
            printf("Copying batch results to host\n");
            fflush(stdout);
            // copy batch results to host
            HANDLE_ERROR(cudaMemcpy(
                &corrs[batch_start * (corr_width + num_phen)],
                gpu_corrs,
                (corr_width + num_phen) * batch_size * sizeof(float),
                cudaMemcpyDeviceToHost
            ));
            batch_start = row_ix + 1;
            fflush(stdout);
        }

        if (row_ix >= last_full_row)
        {
            curr_width -= 1;
            max_row_in += 1;
        }

        if (row_in == (max_row_in - 1))
        {
            printf("At marker #%u;\n", row_ix);
            printf("Loading new marker batch \n");
            fflush(stdout);
            // get new marker data batch
            HANDLE_ERROR(cudaMemcpy(
                gpu_marker_vals,
                &marker_vals[genotype_col_bytes * (row_ix + 1)],
                batch_marker_vals_bytes,
                cudaMemcpyHostToDevice
            ));
            row_in = 0;
        }
        else
        {
            row_in += 1;
        }
    }

    printf("Finished computing marker-marker and marker-phen \n");

    // copy remaining results to host
    if ((num_markers % batch_size) != 0)
    {
        size_t last_batch_size = num_markers - batch_start;
        HANDLE_ERROR(cudaMemcpy(
            &corrs[batch_start * (corr_width + num_phen)],
            gpu_corrs,
            (corr_width + num_phen) * last_batch_size * sizeof(float),
            cudaMemcpyDeviceToHost
        ));
    }

    printf("Computing phen-phen corrs \n");
    fflush(stdout);
    // phen-phen corr
    size_t num_phen_corrs = num_phen * (num_phen - 1) / 2;
    BLOCKS_PER_GRID = dim3(num_phen_corrs, 1, 1);
    THREADS_PER_BLOCK = dim3(NUMTHREADS, 1, 1);
    float *gpu_phen_corrs;
    HANDLE_ERROR(cudaMalloc(&gpu_phen_corrs, num_phen_corrs * sizeof(float)));
    phen_corr_pearson_scan<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
        gpu_phen_vals, num_individuals, num_phen, gpu_phen_corrs
    );
    CudaCheckError();
    std::vector<float> phen_corrs_tmp(num_phen_corrs, 0.0);
    HANDLE_ERROR(cudaMemcpy(
        phen_corrs_tmp.data(),
        gpu_phen_corrs,
        num_phen_corrs * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    size_t num_row_entries = num_phen;
    size_t lin_ix = 0;
    for (size_t pix = 0; pix < num_phen; ++pix)
    {
        for (size_t vix = pix + 1; vix < num_row_entries; ++vix)
        {
            size_t pix1 = pix;
            size_t pix2 = vix;
            corrs[(num_markers + pix1) * corr_row_len + corr_width + pix2] = phen_corrs_tmp[lin_ix];
            corrs[(num_markers + pix2) * corr_row_len + corr_width + pix1] = phen_corrs_tmp[lin_ix];
            lin_ix += 1;
        }
    }

    HANDLE_ERROR(cudaFree(gpu_phen_corrs));
    HANDLE_ERROR(cudaFree(gpu_corrs));
    HANDLE_ERROR(cudaFree(gpu_marker_vals));
    HANDLE_ERROR(cudaFree(gpu_marker_mean));
    HANDLE_ERROR(cudaFree(gpu_marker_std));
    HANDLE_ERROR(cudaFree(gpu_phen_vals));
}
