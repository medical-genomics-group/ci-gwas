#pragma once

#include <mps/bfiles_base.h>
#include <mps/io.h>

#include <string>

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