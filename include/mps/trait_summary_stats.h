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
    int num_phen;

   public:
    TraitSummaryStats(const std::string path);
};