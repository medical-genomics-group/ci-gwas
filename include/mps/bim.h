#pragma once
#include <string>
#include <vector>

#define BIM_NUM_COLS 6

struct bimInfo
{
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
};

void split_bim_line(std::string line, std::string *buf);

auto prep_bim(std::string bim_path, std::string out_dir) -> bimInfo;

int median(std::vector<int> &v);

auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t;