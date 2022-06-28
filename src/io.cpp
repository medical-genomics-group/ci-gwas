// Generic functions for reading and writing 

#include <mps/io.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>

auto count_lines(std::string file_path) -> int {
     int number_of_lines = 0;
     std::string line;
     std::ifstream fin(file_path);
 
     while (std::getline(fin, line))
         ++number_of_lines;
 
     return number_of_lines;
}

auto make_path(const std::string out_dir, const std::string chr_id, const std::string suffix) -> std::string {
    std::string filename = chr_id + suffix;
    std::string outpath;
    if (out_dir.back() != '/') {
        outpath = out_dir + '/' + filename;
    } else {
        outpath = out_dir + filename;
    }
    return outpath;
}

// TODO: I don't ever need dest to have dynamic size, this should accept arrays
void read_n_bytes_from_binary(const std::string path, const size_t nbytes, std::vector<unsigned char> &dest)
{
    std::ifstream bin_file(path, std::ios::binary);
    bin_file.read(reinterpret_cast<char*>(dest.data()), nbytes); 
}

void read_floats_from_lines(const std::string path, std::vector<float> &dest)
{
    std::string line;
    std::ifstream fin(path);
    
    while (std::getline(fin, line)) {
        dest.push_back(std::stof(line));
    }
}

auto read_floats_from_lines(const std::string path) -> std::vector<float>
{
    std::string line;
    std::ifstream fin(path);
    std::vector<float> res = {};

    while (std::getline(fin, line)) {
        res.push_back(std::stof(line));
    }

    return res;
}

auto read_ints_from_lines(const std::string path) -> std::vector<int>
{
    std::string line;
    std::ifstream fin(path);
    std::vector<int> res = {};

    while (std::getline(fin, line)) {
        res.push_back(std::stoi(line));
    }

    return res;
}


void write_bed(
        const std::vector<unsigned char> &out_buf,
        const std::string out_dir,
        const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".bed");
    std::ofstream bedout;
    bedout.open(outpath, std::ios::out|std::ios::binary);
    std::copy(out_buf.cbegin(), out_buf.cend(), std::ostream_iterator<unsigned char>(bedout));
    bedout.close();
}

void write_means(
        const std::vector<float> &chr_marker_means,
        const std::string out_dir,
        const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".means");
    std::ofstream fout;
    fout.open(outpath, std::ios::out);

    for (size_t i = 0; i < chr_marker_means.size(); ++i) {
        fout << chr_marker_means[i] << std::endl;
    }

    fout.close();
}

void write_stds(
        const std::vector<float> &chr_marker_stds,
        const std::string out_dir,
        const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".stds");
    std::ofstream fout;
    fout.open(outpath, std::ios::out);

    for (size_t i = 0; i < chr_marker_stds.size(); ++i) {
        fout << chr_marker_stds[i] << std::endl;
    }

    fout.close();
}

void write_dims(
        const size_t num_individuals,
        const size_t num_markers,
        const std::string out_dir,
        const std::string chr_id)
{
    std::string outpath = make_path(out_dir, chr_id, ".dims");
    std::ofstream fout;
    fout.open(outpath, std::ios::out);
    fout << num_individuals << std::endl;
    fout << num_markers << std::endl;
    fout.close();
}

