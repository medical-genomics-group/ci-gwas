#pragma once

void cu_marker_corr_pearson_npn(const unsigned char *marker_vals, const size_t num_markers,
                                const size_t num_individuals, float *marker_corrs);

void cu_marker_corr_pearson_npn_batched(const unsigned char *marker_vals, const size_t num_markers,
                                        const size_t num_individuals, const size_t row_set_size,
                                        float *marker_corrs);

void cu_marker_corr_pearson(const unsigned char *marker_vals, const size_t num_markers,
                            const size_t num_individuals, const float *marker_mean,
                            const float *marker_std, float *marker_corrs);

void cu_marker_corr_pearson_batched(const unsigned char *marker_vals, const size_t num_markers,
                                    const size_t num_individuals, const float *marker_mean,
                                    const float *marker_std, const size_t batch_stripe_width,
                                    float *marker_corrs);

void cu_corr_pearson_npn(const unsigned char *marker_vals, const float *phen_vals,
                         const size_t num_markers, const size_t num_individuals,
                         const size_t num_phen, const float *marker_mean, const float *marker_std,
                         float *marker_corrs, float *marker_phen_corrs, float *phen_corrs);

void cu_corr_pearson_npn_batched(const unsigned char *marker_vals, const float *phen_vals,
                                 const size_t num_markers, const size_t num_individuals,
                                 const size_t num_phen, const float *marker_mean,
                                 const float *marker_std, const size_t row_width,
                                 float *marker_corrs, float *marker_phen_corrs, float *phen_corrs);