#pragma once

#define NUMTHREADS 16

void cu_corr_npn(const unsigned char *a, const size_t num_markers, const size_t num_individuals,
                 float *results);

__global__ void cu_marker_corr_npn(const unsigned char *a, const size_t num_markers,
                                   const size_t num_individuals, float *results);