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
};