#pragma once

#include <unordered_set>
#include <vector>

std::unordered_set<int> parent_set(
    const std::vector<int> &G, const size_t num_var, const size_t num_markers, const int max_depth
);

std::vector<int> set_to_vec(const std::unordered_set<int>);
