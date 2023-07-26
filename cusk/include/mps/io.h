#pragma once

#include <mps/bfiles_base.h>
#include <mps/bim.h>
#include <mps/marker_block.h>

#include <string>
#include <vector>

float zero_if_nan(float value);

void split_line(std::string line, std::string *buf, size_t ncols);

std::vector<std::string> split_line(std::string line);

int count_lines(const std::string file_path);

class BedDims
{
   private:
    size_t num_samples;
    size_t num_markers;

   public:
    BedDims(size_t num_samples, size_t num_markers)
        : num_samples(num_samples), num_markers(num_markers)
    {
    }

    BedDims(std::string path)
    {
        std::string dim_line[2];
        std::string line;
        std::ifstream dim_file(path);
        std::getline(dim_file, line);
        split_line(line, dim_line, 2);
        num_samples = std::stoi(dim_line[0]);
        num_markers = std::stoi(dim_line[1]);
    }

    BedDims(BfilesBase bfiles)
    {
        num_samples = count_lines(bfiles.fam());
        num_markers = count_lines(bfiles.bim());
    }

    bool operator==(const BedDims &other) const
    {
        return (this->num_samples == other.num_samples) && (this->num_markers == other.num_markers);
    }

    size_t get_num_markers() const { return num_markers; }

    size_t get_num_samples() const { return num_samples; }

    size_t bytes_per_col() const { return (num_samples + 3) / 4; }

    void to_file(std::string path) const
    {
        std::ofstream fout;
        fout.open(path, std::ios::out);
        fout << num_samples << "\t" << num_markers << std::endl;
        fout.close();
    }
};

std::vector<float> read_floats_from_line_range(const std::string path, size_t first, size_t last);

std::vector<unsigned char> read_block_from_bed(
    std::string path, MarkerBlock block, BedDims dims, BimInfo bim
);

std::vector<unsigned char> read_chr_from_bed(
    std::string path, std::string chr_id, BimInfo bim, BedDims dim
);

std::string make_path(std::string out_dir, std::string chr_id, std::string suffix);

std::vector<int> read_ints_from_lines(std::string path);

void read_floats_from_lines(std::string path, std::vector<float> &dest);

std::vector<float> read_floats_from_lines(std::string path);

void write_marker_blocks_to_file(const std::vector<MarkerBlock> &blocks, const std::string path);

void write_bed(
    const std::vector<unsigned char> &out_buf, const std::string out_dir, const std::string chr_id
);

void write_means(
    const std::vector<float> &chr_marker_means, const std::string out_dir, const std::string chr_id
);

void write_stds(
    const std::vector<float> &chr_marker_stds, const std::string out_dir, const std::string chr_id
);

void write_single_column_file_with_suffix(
    const std::vector<float> &data, const std::string outpath
);

void write_single_column_file_with_suffix(const std::vector<int> &data, const std::string outpath);

void write_single_column_file_with_suffix(
    const std::vector<float> &data,
    const std::string out_dir,
    const std::string file_stem,
    const std::string suffix
);

void write_dims(
    const size_t num_individuals,
    const size_t num_markers,
    const std::string out_dir,
    const std::string chr_id
);

void read_n_bytes_from_binary(
    const std::string path, const size_t nbytes, std::vector<unsigned char> &dest
);

std::vector<float> read_floats_from_binary(const std::string path);

void write_ints_to_binary(const int *data, const size_t nvals, const std::string path);

void write_floats_to_binary(const float *data, const size_t nvals, const std::string path);

void append_floats_to_binary(const float *data, const size_t nvals, const std::string path);

std::vector<MarkerBlock> read_blocks_from_file(const std::string path);

void check_prepped_bed_path(const std::string basepath);

void check_bed_path(const std::string basepath);

void check_path(const std::string path);

bool path_exists(const std::string path);

std::vector<float> read_correlations_from_mtx(const std::string path);