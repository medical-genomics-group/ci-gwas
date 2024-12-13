#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class TraitSummaryStats
{
   private:
    std::vector<std::string> header;
    std::vector<float> corrs;
    std::vector<float> sample_sizes;
    int num_phen;

   public:
    TraitSummaryStats(const std::string corr_path);
    TraitSummaryStats(const std::string corr_path, const std::string se_path);
    TraitSummaryStats(const std::string corr_path, const float sample_size);
    std::vector<std::string> get_header() { return header; };
    std::vector<float> get_sample_sizes() { return sample_sizes; };
    std::vector<float> get_corrs() { return corrs; };
    int get_num_phen() { return num_phen; };
};