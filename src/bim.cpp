#include <mps/bim.h>
#include <mps/io.h>
#include <mps/prep.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void split_bim_line(std::string line, std::string *buf) { split_line(line, buf, BIM_NUM_COLS); }

int median(std::vector<int> &v)
{
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

// TODO: write function that deteremines number of markers within a given distance
// this assumes that the bim contains only a single chromosome
auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t
{
    std::string bim_line[6];
    std::string line;
    std::ifstream bim(bim_path);

    std::vector<unsigned long> positions;
    std::vector<int> marker_nums;

    size_t pa = 0;
    size_t pb = 0;
    while (std::getline(bim, line))
    {
        split_bim_line(line, bim_line);
        positions.push_back(std::stoul(bim_line[3]));
        while ((positions[pb] - positions[pa]) > dthr)
        {
            marker_nums.push_back(pb - pa - 1);
            pa += 1;
        }
        pb += 1;
    }

    return (size_t)median(marker_nums);
}