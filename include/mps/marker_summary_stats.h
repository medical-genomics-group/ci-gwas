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
    MarkerSummaryStats(const std::string path);
    std::vector<float> get_corrs() { return corrs; };
    int get_num_markers() { return num_markers; };
};