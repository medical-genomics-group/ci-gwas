#pragma once

#include <vector>

std::vector<int> parent_set(
    const std::vector<int> &G,
    const size_t num_var,
    const size_t num_markers,
    const int max_depth);
