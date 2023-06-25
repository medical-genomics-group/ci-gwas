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
    int num_phen;
    int num_markers;

   public:
    MarkerTraitSummaryStats(const std::string path);
};