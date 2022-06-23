#pragma once
#include <string>

struct bimInfo {
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
};

void split_bim_line(std::string line, std::string *buf);
auto count_lines(std::string file_path) -> int;
auto parse_bim(std::string bim_path) -> bimInfo;

