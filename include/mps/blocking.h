#pragma once

#include <mps/io.h>

#include <vector>

std::vector<MarkerBlock> block_chr(
    const std::vector<float> &forward_corr_sums, const std::string chr_id, const int max_block_size
);

std::vector<double> hanning_smoothing(const std::vector<float> &v, const int window_size);
