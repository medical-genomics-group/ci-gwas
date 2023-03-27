#pragma once
#include <mps/marker_block.h>

#include <cmath>
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

    BimInfo(std::string path);

    void check_chr_id(const std::string chr_id) const;

    size_t get_num_markers_on_chr(const std::string chr_id)
    {
        check_chr_id(chr_id);
        return num_markers_on_chr[chr_id2ix[chr_id]];
    }

    size_t get_global_chr_start(const std::string chr_id)
    {
        check_chr_id(chr_id);
        return global_chr_start[chr_id2ix[chr_id]];
    }

    // inclusive end
    size_t get_global_chr_end(const std::string chr_id)
    {
        check_chr_id(chr_id);
        size_t chr_ix = chr_id2ix[chr_id];
        return global_chr_start[chr_ix] + num_markers_on_chr[chr_ix] - 1;
    }

    std::vector<MarkerBlock> equisized_blocks(const int size) const
    {
        std::vector<MarkerBlock> res;
        for (size_t i = 0; i < chr_ids.size(); i++)
        {
            size_t num_blocks(num_markers_on_chr[i] + size - 1) / size;
            for (size_t bid = 0; bid < num_blocks; bid++)
            {
                size_t bstart = bid * size;
                size_t bend = min(bstart + size - 1, num_markers_on_chr[i]);
                res.push_back(MarkerBlock(chr_ids[i], bstart, bend));
            }
        }
        return res;
    }
};

int median(std::vector<int> &v);

auto num_markers_within_distance(std::string bim_path, unsigned long dthr) -> size_t;
