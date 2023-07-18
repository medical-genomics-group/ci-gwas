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
    std::vector<std::string> get_header() { return header; };
    std::vector<float> get_corrs() { return corrs; };
    int get_num_phen() { return num_phen; };
};