#pragma once
#include <string>
#include <vector>

struct bimInfo
{
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
};

void prep_bed(
    std::string bed_path,
    std::string bim_path,
    std::string fam_path,
    std::string out_dir,
    size_t mem_gb
);

void split_bim_line(std::string line, std::string *buf);

auto parse_bim(std::string bim_path) -> bimInfo;
