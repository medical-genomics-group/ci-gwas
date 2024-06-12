#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class MarkerTraitSummaryStats
{
   private:
    std::vector<std::string> header;
    std::vector<float> corrs;
    std::vector<float> sample_sizes;
    int num_phen;
    int num_markers;

   public:
    MarkerTraitSummaryStats(const std::string path);
    MarkerTraitSummaryStats(const std::string path, const MarkerBlock block);
    MarkerTraitSummaryStats(const std::string corr_path, const std::string se_path, const MarkerBlock block, const float nan_sample_size_fill_val);
    MarkerTraitSummaryStats(const std::string path, const std::vector<int> marker_ixs);
    std::vector<std::string> get_header() { return header; };
    std::vector<float> get_corrs() { return corrs; };
    std::vector<float> get_sample_sizes() { return sample_sizes; };
    int get_num_phen() { return num_phen; };
    int get_num_markers() { return num_markers; };
};