#pragma once

#include <string>
#include <vector>

const int BED_PREFIX_BYTES = 3;

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

    size_t get_num_markers() const
    {
        return num_markers;
    }

    size_t get_num_samples() const
    {
        return num_samples;
    }

    size_t bytes_per_col() const
    {
        return (num_samples + 3) / 4;
    }
};

class MarkerBlock
{
  private:
    size_t first_marker_ix;
    size_t last_marker_ix;

  public:
    MarkerBlock(size_t ix1, size_t ix2)
        : first_marker_ix(ix1), last_marker_ix(ix2)
    {
    }

    bool operator==(const MarkerBlock &other) const
    {
        return (this->first_marker_ix == other.first_marker_ix) && (this->last_marker_ix == other.last_marker_ix);
    }

    size_t get_first_marker_ix() const
    {
        return first_marker_ix;
    }

    size_t get_last_marker_ix() const
    {
        return last_marker_ix;
    }

    size_t block_size() const
    {
        return last_marker_ix - first_marker_ix + 1;
    }
};

std::vector<unsigned char> read_block_from_bed(
    std::string path,
    MarkerBlock block,
    BedDims dims);

void split_line(
    std::string line,
    std::string *buf,
    size_t ncols);

std::vector<std::string> split_line(
    std::string line);

std::string make_path(
    std::string out_dir,
    std::string chr_id,
    std::string suffix);

std::vector<int> read_ints_from_lines(
    std::string path);

void read_floats_from_lines(
    std::string path,
    std::vector<float> &dest);

std::vector<float> read_floats_from_lines(
    std::string path);

void write_bed(
    const std::vector<unsigned char> &out_buf,
    const std::string out_dir,
    const std::string chr_id);

void write_means(
    const std::vector<float> &chr_marker_means,
    const std::string out_dir,
    const std::string chr_id);

void write_stds(
    const std::vector<float> &chr_marker_stds,
    const std::string out_dir,
    const std::string chr_id);

void write_single_column_file_with_suffix(
    const std::vector<float> &data,
    const std::string out_dir,
    const std::string file_stem,
    const std::string suffix);

int count_lines(
    const std::string file_path);

void write_dims(
    const size_t num_individuals,
    const size_t num_markers,
    const std::string out_dir,
    const std::string chr_id);

void read_n_bytes_from_binary(
    const std::string path,
    const size_t nbytes,
    std::vector<unsigned char> &dest);

std::vector<float> read_floats_from_binary(
    const std::string path);

void write_floats_to_binary(
    const float *data,
    const size_t nvals,
    const std::string path);

std::vector<MarkerBlock> read_blocks_from_file(
    const std::string path);