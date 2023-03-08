#pragma once
#include <string>
#include <vector>

struct Phen
{
    size_t num_samples;
    size_t num_phen;
    std::vector<float> data;

    size_t get_num_samples() const
    {
        return num_samples;
    }

    size_t get_num_phen() const
    {
        return num_phen;
    }
};

Phen load_phen(std::string phen_path);