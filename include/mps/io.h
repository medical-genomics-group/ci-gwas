#pragma once

#include <string>
#include <vector>

void split_line(std::string line, std::string *buf, size_t ncols);

auto split_line(std::string line) -> std::vector<std::string>;

auto make_path(std::string out_dir, std::string chr_id, std::string suffix) -> std::string;

auto read_ints_from_lines(std::string path) -> std::vector<int>;

void read_floats_from_lines(std::string path, std::vector<float> &dest);

auto read_floats_from_lines(std::string path) -> std::vector<float>;

void write_bed(
    const std::vector<unsigned char> &out_buf, const std::string out_dir, const std::string chr_id);

void write_means(
    const std::vector<float> &chr_marker_means, const std::string out_dir, const std::string chr_id);

void write_stds(
    const std::vector<float> &chr_marker_stds, const std::string out_dir, const std::string chr_id);

auto count_lines(std::string file_path) -> int;

void write_dims(
    const size_t num_individuals,
    const size_t num_markers,
    const std::string out_dir,
    const std::string chr_id);

void read_n_bytes_from_binary(
    const std::string path, const size_t nbytes, std::vector<unsigned char> &dest);

auto read_floats_from_binary(const std::string path) -> std::vector<float>;

void write_floats_to_binary(const float *data, const size_t nvals, const std::string path);
