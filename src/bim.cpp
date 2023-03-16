#include <mps/bim.h>
#include <mps/io.h>
#include <mps/prep.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

BimInfo::BimInfo(std::string path)
{
    number_of_lines = 0;
    chr_ids = {};
    num_markers_on_chr = {};
    chr_id2ix = {};
    global_chr_start = {};

    std::string bim_line[BIM_NUM_COLS];
    std::string line;
    std::ifstream bim(path);

    size_t chr_ix = 0;

    while (std::getline(bim, line))
    {
        split_bim_line(line, bim_line);
        if ((number_of_lines == 0) || (bim_line[0] != chr_ids.back()))
        {
            global_chr_start.push_back(number_of_lines);
            chr_id2ix[bim_line[0]] = chr_ix;
            chr_ids.push_back(bim_line[0]);
            num_markers_on_chr.push_back(0);
            ++chr_ix;
        }
        ++num_markers_on_chr.back();
        ++number_of_lines;
    }
}

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