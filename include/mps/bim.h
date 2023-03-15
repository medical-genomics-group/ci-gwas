#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#define BIM_NUM_COLS 6

struct BimInfo
{
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
    std::unordered_map<std::string, size_t> chr_id2ix;

    BimInfo(std::string path)
    {
        number_of_lines = 0;
        chr_ids = {};
        num_markers_on_chr = {};
        chr_id2ix = {};

        std::string bim_line[BIM_NUM_COLS];
        std::string line;
        std::ifstream bim(path);

        size_t chr_ix = 0;

        while (std::getline(bim, line))
        {
            split_bim_line(line, bim_line);
            if ((number_of_lines == 0) || (bim_line[0] != chr_ids.back()))
            {
                chr_id2ix[bim_line[0]] = chr_ix;
                chr_ids.push_back(bim_line[0]);
                num_markers_on_chr.push_back(0);
                ++chr_ix;
            }
            ++num_markers_on_chr.back();
            ++number_of_lines;
        }
    }

    size_t num_markers_on_chr(std::string chr_id) const
    {
        if (!chr_id2ix.contains(chr_id))
        {
            std::cerr << "chr id does not match .bim content" << std::endl;
            exit(1);
        }
        return chr_id2ix[chr_id];
    }

    size_t chr_start_global_ix(std::string chr_id) const
    {
        size_t res = 0;
        for (size_t i = 0; i < chr_ids.size(); i++)
        {
            if (chr_id == chr_ids[i])
            {
                return res;
            }
            res += num_markers_on_chr[i];
        }
        std::cerr << "chr id does not match .bim content" << std::endl;
        exit(1);
    }
};

void split_bim_line(std::string line, std::string *buf);

int median(std::vector<int> &v);

auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t;
