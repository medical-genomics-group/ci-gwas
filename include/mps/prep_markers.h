#pragma once
#include <string>

struct bimInfo {
    size_t number_of_lines;
    std::vector<std::string> chr_ids;
    std::vector<size_t> num_markers_on_chr;
};


void write_bed(unsigned char *out_buf, std::string out_dir, std::string chr_id);
void write_means(float *chr_marker_means, size_t num_vals, std::string out_dir, std::string chr_id);
void write_stds(float *chr_marker_stds, size_t num_vals, std::string out_dir, std::string chr_id);
void prep_bed(std::string bed_path, std::string bim_path, std::string fam_path, std::string out_dir, size_t mem_gb);
void split_bim_line(std::string line, std::string *buf);
auto count_lines(std::string file_path) -> int;
auto parse_bim(std::string bim_path) -> bimInfo;

