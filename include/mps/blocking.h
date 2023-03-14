#pragma once

#include <mps/io.h>
#include <vector>

std::vector<MarkerBlock> block_chr(
    std::vector<float> &diag_sums,
    std::string chr_id,
    int max_block_size);