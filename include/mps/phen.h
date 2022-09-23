#pragma once
#include <string>
#include <vector>

struct Phen
{
    size_t num_samples;
    size_t num_phen;
    std::vector<float> data;
};

auto load_phen(std::string phen_path) -> Phen;