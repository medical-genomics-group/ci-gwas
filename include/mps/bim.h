#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define BIM_NUM_COLS 6

struct bimInfo
{
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;

    bimInfo(std::string path)
    {
        number_of_lines = 0;
        chr_ids = {};
        num_markers_on_chr = {};

        std::string bim_line[BIM_NUM_COLS];
        std::string line;
        std::ifstream bim(path);

        while (std::getline(bim, line))
        {
            split_bim_line(line, bim_line);
            if ((number_of_lines == 0) || (bim_line[0] != chr_ids.back()))
            {
                chr_ids.push_back(bim_line[0]);
                num_markers_on_chr.push_back(0);
            }
            ++num_markers_on_chr.back();
            ++number_of_lines;
        }
    }
};

void split_bim_line(std::string line, std::string *buf);

auto prep_bim(std::string bim_path, std::string out_dir) -> bimInfo;

int median(std::vector<int> &v);

auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t;
