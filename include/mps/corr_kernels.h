#pragma once

#define NUMTHREADS 64

__global__ void ix_from_linear(
    const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix
);

__global__ void bed_marker_corr_pearson_batched(
    const unsigned char *row_marker_vals,
    const unsigned char *col_marker_vals,
    const size_t num_col_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean_row,
    const float *marker_std_row,
    const float *marker_mean_col,
    const float *marker_std_col,
    float *results
);

__global__ void bed_marker_corr_pearson_batched_row(
    const unsigned char *row_marker_vals,
    const size_t num_rows,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
);

__global__ void bed_marker_phen_corr_pearson(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
);

__global__ void phen_corr_pearson(
    const float *phen_vals, const size_t num_individuals, const size_t num_phen, float *results
);

__global__ void bed_marker_corr_kendall_npn(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
);

__global__ void bed_marker_corr_pearson(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
);

__global__ void bed_marker_corr_kendall_npn_batched(
    const unsigned char *row_marker_vals,
    const unsigned char *col_marker_vals,
    const size_t num_col_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
);

__global__ void bed_marker_corr_kendall_npn_batched_row(
    const unsigned char *row_marker_vals,
    const size_t num_rows,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
);