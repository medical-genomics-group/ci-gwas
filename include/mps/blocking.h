#pragma once

#include <mps/io.h>
#include <vector>

std::vector<MarkerBlock> block_chr(
    const std::vector<float> &diag_sums,
    const std::string chr_id,
    const int max_block_size);
