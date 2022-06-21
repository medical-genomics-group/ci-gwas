#pragma once

#define NUMTHREADS 16

void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results);

__global__ void cu_marker_corr_npn(const unsigned char *a, const size_t num_markers,
                                   const size_t num_individuals, float *results);

void cu_corr_npn(const unsigned char *marker_vals, const float *phen_vals, const size_t num_markers,
                 const size_t num_individuals, const size_t num_phen, const float marker_mean,
                 const float marker_std, float *results);

__global__ void bed_marker_corr_npn(const unsigned char *a, const size_t num_markers,
                                    const size_t num_individuals, const size_t col_len_bytes,
                                    float *results);
