#pragma once

#include <mps/io.h>

#include <vector>

std::vector<float> cal_mcorrk(
    const std::vector<unsigned char> &bed_vals,
    const BedDims dim,
    const size_t num_markers,
    const size_t device_mem_gb
);

std::vector<float> cal_mcorrk_banded(
    const std::vector<unsigned char> &bed_vals,
    const BedDims dim,
    const size_t num_markers,
    const size_t corr_width,
    const float device_mem_gb
);

std::vector<float> marker_corr_banded_mat_row_abs_sums(
    const std::vector<float> &marker_corrs, size_t corr_width
);

std::vector<float> marker_corr_mat_antidiag_sums(const std::vector<float> &marker_corrs);

void marker_corr_mat_antidiag_sums(
    const size_t num_markers, const float *marker_corrs, float *sums
);

void cu_ix_from_linear(const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix);

void cu_phen_corr_pearson_npn(
    const float *phen_vals, const size_t num_individuals, const size_t num_phen, float *phen_corrs
);

void cu_marker_phen_corr_pearson(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    float *marker_phen_corrs
);

void cu_marker_corr_pearson_npn(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    float *marker_corrs
);

void cu_marker_corr_pearson_npn_batched(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t row_set_size,
    float *marker_corrs
);

void cu_marker_corr_pearson_npn_batched_sparse(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t corr_width,
    const size_t batch_size,
    float *corrs
);

void cu_marker_corr_pearson(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const float *marker_mean,
    const float *marker_std,
    float *marker_corrs
);

void cu_marker_corr_pearson_batched(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const float *marker_mean,
    const float *marker_std,
    const size_t batch_stripe_width,
    float *marker_corrs
);

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
);

void cu_corr_pearson_npn_batched(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const float *marker_mean,
    const float *marker_std,
    const size_t row_width,
    float *marker_corrs,
    float *marker_phen_corrs,
    float *phen_corrs
);

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
);
