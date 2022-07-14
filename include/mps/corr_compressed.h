#pragma once

#define NUMTHREADS 128

void cu_corr_npn(const unsigned char *marker_vals, const float *phen_vals, const size_t num_markers,
                 const size_t num_individuals, const size_t num_phen, const float *marker_mean,
                 const float *marker_std, float *marker_corrs, float *marker_phen_corrs,
                 float *phen_corrs);

__global__ void marker_phen_corr_pearson(const unsigned char *marker_vals, const float *phen_vals,
                                         const size_t num_markers, const size_t num_individuals,
                                         const size_t num_phen, const size_t col_len_bytes,
                                         const float *marker_mean, const float *marker_std,
                                         float *results);

__global__ void phen_corr_pearson(const float *phen_vals, const size_t num_individuals,
                                  const size_t num_phen, float *results);

__global__ void bed_marker_corr_kendall_npn(const unsigned char *marker_vals,
                                            const size_t num_markers, const size_t num_individuals,
                                            const size_t col_len_bytes, float *results);

__global__ void marker_corr_pearson(const unsigned char *marker_vals,
                                     const size_t num_markers,
                                     const size_t num_individuals,
                                     const size_t col_len_bytes,
                                     const float *marker_mean,
                                     const float *marker_std,
                                     float *results);
