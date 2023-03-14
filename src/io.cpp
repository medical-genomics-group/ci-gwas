// Generic functions for reading and writing

#include <mps/io.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <sys/stat.h>

std::vector<std::string> split_line(std::string line)
{
    std::vector<std::string> res;
    std::istringstream ss(line);
    std::string word;
    while (ss >> word)
    {
        res.push_back(word);
    }
    return res;
}

void split_line(std::string line, std::string *buf, size_t ncols)
{
    std::istringstream ss(line);

    std::string word;
    size_t word_ix = 0;

    while (ss >> word)
    {
        assert((word_ix < ncols) && "too many entries in line, aborting.");
        buf[word_ix] = word;
        word_ix += 1;
    }
}

int count_lines(const std::string file_path)
{
    int number_of_lines = 0;
    std::string line;
    std::ifstream fin(file_path);

    while (std::getline(fin, line))
        ++number_of_lines;

    return number_of_lines;
}

auto make_path(const std::string out_dir, const std::string file_stem, const std::string suffix)
    -> std::string
{
    std::string filename = file_stem + suffix;
    std::string outpath;
    if (out_dir.empty())
    {
        outpath = out_dir + filename;
    }
    else if (out_dir.back() != '/')
    {
        outpath = out_dir + '/' + filename;
    }
    else
    {
        outpath = out_dir + filename;
    }
    return outpath;
}

std::vector<MarkerBlock> read_blocks_from_file(const std::string path)
{
    std::string block_line[3];
    std::string line;
    std::ifstream block_file(path);
    std::vector<MarkerBlock> blocks;

    while (std::getline(block_file, line))
    {
        split_line(line, block_line, 3);
        blocks.push_back(MarkerBlock((std::string)block_line[0], std::stoi(block_line[1]), std::stoi(block_line[2])));
    }

    return blocks;
}

// TODO: I don't ever need dest to have dynamic size, this should accept arrays
void read_n_bytes_from_binary(
    const std::string path, const size_t nbytes, std::vector<unsigned char> &dest)
{
    dest.resize(nbytes);
    std::ifstream bin_file(path, std::ios::binary);
    bin_file.read(reinterpret_cast<char *>(dest.data()), nbytes);
}

void read_floats_from_lines(const std::string path, std::vector<float> &dest)
{
    std::string line;
    std::ifstream fin(path);

    while (std::getline(fin, line))
    {
        dest.push_back(std::stof(line));
    }
}

std::vector<float> read_floats_from_lines(const std::string path)
{
    std::string line;
    std::ifstream fin(path);
    std::vector<float> res = {};

    while (std::getline(fin, line))
    {
        res.push_back(std::stof(line));
    }

    return res;
}

std::vector<float> read_floats_from_line_range(const std::string path, size_t first, size_t last)
{
    std::string line;
    std::ifstream fin(path);
    std::vector<float> res = {};

    size_t line_ix = 0;
    while (std::getline(fin, line))
    {
        if (line_ix > last)
        {
            break;
        }
        else if (line_ix >= first)
        {
            res.push_back(std::stof(line));
        }
        ++line_ix;
    }

    return res;
}

std::vector<int> read_ints_from_lines(const std::string path)
{
    std::string line;
    std::ifstream fin(path);
    std::vector<int> res = {};

    while (std::getline(fin, line))
    {
        res.push_back(std::stoi(line));
    }

    return res;
}

std::vector<float> read_floats_from_binary(const std::string path)
{
    float f;
    std::ifstream fin(path, std::ios::binary);
    std::vector<float> res = {};

    while (fin.read(reinterpret_cast<char *>(&f), sizeof(float)))
        res.push_back(f);

    return res;
}

std::vector<unsigned char> read_block_from_bed(
    std::string path,
    MarkerBlock block,
    BedDims dims,
    BimInfo bim)
{
    size_t chr_start = bim.chr_start_global_ix(block.get_chr_id());
    size_t block_bytes = dims.bytes_per_col() * block.block_size();
    std::vector<unsigned char> res(block_bytes, 0);
    std::ifstream fin(path, std::ios::binary);
    fin.seekg(BED_PREFIX_BYTES + chr_start + dims.bytes_per_col() * block.get_first_marker_ix());
    fin.read(reinterpret_cast<char *>(res.data()), block_bytes);
    return res;
}

void write_floats_to_binary(const float *data, const size_t nvals, const std::string path)
{
    std::ofstream fout;
    fout.open(path, std::ios::out | std::ios::binary);
    for (size_t i = 0; i < nvals; ++i)
    {
        fout.write(reinterpret_cast<const char *>(&data[i]), sizeof(float));
    }
    fout.close();
}

void write_bed(
    const std::vector<unsigned char> &out_buf, const std::string out_dir, const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".bed");
    std::ofstream bedout;
    bedout.open(outpath, std::ios::out | std::ios::binary);
    std::copy(out_buf.cbegin(), out_buf.cend(), std::ostream_iterator<unsigned char>(bedout));
    bedout.close();
}

void write_means(
    const std::vector<float> &chr_marker_means, const std::string out_dir, const std::string chr_id)
{
    write_single_column_file_with_suffix(chr_marker_means, out_dir, chr_id, ".means");
}

void write_stds(
    const std::vector<float> &chr_marker_stds, const std::string out_dir, const std::string chr_id)
{
    write_single_column_file_with_suffix(chr_marker_stds, out_dir, chr_id, ".stds");
}

void write_single_column_file_with_suffix(
    const std::vector<float> &data, const std::string out_dir, const std::string file_stem, const std::string suffix)
{
    std::string outpath = make_path(out_dir, file_stem, suffix);
    std::ofstream fout;
    fout.open(outpath, std::ios::out);

    std::cout << "Writing to: " << outpath << std::endl;

    for (size_t i = 0; i < data.size(); ++i)
    {
        fout << data[i] << std::endl;
    }

    fout.close();
}

void write_dims(
    const size_t num_individuals,
    const size_t num_markers,
    const std::string out_dir,
    const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".dim");
    std::ofstream fout;
    fout.open(outpath, std::ios::out);
    fout << num_individuals << std::endl;
    fout << num_markers << std::endl;
    fout.close();
}

void check_path(const std::string path)
{
    if (!path_exists(path))
    {
        std::cout << "file or directory not found: " << path << std::endl;
        exit(1);
    }
}

void check_prepped_bed_path(const std::string basepath)
{
    std::string req_suffixes[6] = {".bed", ".dim", ".means", ".stds", ".bim", ".fam"};
    for (auto &suffix : req_suffixes)
    {
        check_path(basepath + suffix);
    }
}

void check_bed_path(const std::string basepath)
{
    std::string req_suffixes[6] = {".bed", ".bim", ".fam"};
    for (auto &suffix : req_suffixes)
    {
        check_path(basepath + suffix);
    }
}

/**
 * @brief Check if path exists.
 *
 * @param path
 * @return true if path exists
 * @return false if it doesn't
 */
bool path_exists(const std::string path)
{
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}
