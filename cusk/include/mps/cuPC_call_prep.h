#pragma once

#include <vector>

float std_normal_qnorm(const float p);

std::vector<float> threshold_array(const int n, const float alpha);

// Returns the alpha / 2 percentile of the standard normal distribution.
// This differs from `threshold array` in that the percentile is not divided by
// sqrt(n - l - 3). This is because in hetcor, each correlation may have its own
// effective sample size. The full thresholds have to be computed on the fly.
float hetcor_threshold(const float alpha);