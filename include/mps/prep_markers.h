#pragma once
#include <string>
#include <vector>

struct bimInfo {
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
};

auto read_ints_from_lines(std::string path) -> std::vector<int>;

auto read_floats_from_lines(std::string path) -> std::vector<float>;

void write_bed(const std::vector<unsigned char> &out_buf, const std::string out_dir, const std::string chr_id);

void write_means(const std::vector<float> &chr_marker_means, const std::string out_dir, const std::string chr_id);

void write_stds(const std::vector<float> &chr_marker_stds,  const std::string out_dir, const std::string chr_id);

void prep_bed(std::string bed_path, std::string bim_path, std::string fam_path, std::string out_dir, size_t mem_gb);

void split_bim_line(std::string line, std::string *buf);

auto count_lines(std::string file_path) -> int;

auto parse_bim(std::string bim_path) -> bimInfo;

