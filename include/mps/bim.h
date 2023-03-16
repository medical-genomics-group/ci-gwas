#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#define BIM_NUM_COLS 6

void split_bim_line(std::string line, std::string *buf);

struct BimInfo
{
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
    std::vector<size_t> global_chr_start;
    std::unordered_map<std::string, size_t> chr_id2ix;

    BimInfo(std::string path)
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

    size_t check_chr_id(std::string chr_id) const
    {
        if (!chr_id2ix.contains(chr_id))
        {
            std::cerr << "chr id does not match .bim content" << std::endl;
            exit(1);
        }
    }

    size_t get_num_markers_on_chr(std::string chr_id) const
    {
        check_chr_id(chr_id);
        return num_markers_on_chr[chr_id2ix[chr_id]];
    }

    size_t get_global_chr_start(std::string chr_id) const
    {
        check_chr_id(chr_id);
        return global_chr_start[chr_id2ix[chr_id]];
    }

    // inclusive end
    size_t get_global_chr_end(std::chr_id) const
    {
        check_chr_id(chr_id);
        size_t chr_ix = chr_id2ix[chr_id];
        return global_chr_start[chr_ix] + num_markers_on_chr[chr_ix] - 1;
    }
};

int median(std::vector<int> &v);

auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t;
