#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class MarkerSummaryStats
{
   private:
    std::vector<float> corrs;
    int num_markers;

   public:
    MarkerSummaryStats() {}

    MarkerSummaryStats(const std::string path);
    std::vector<float> get_corrs() const { return corrs; };
    int get_num_markers() const { return num_markers; };
};