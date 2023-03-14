#include <mps/bed_lut.h>
#include <mps/bim.h>
#include <mps/io.h>
#include <mps/prep.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

void compute_bed_col_stats_no_impute(
    const std::vector<unsigned char> &bedcol,
    const size_t num_individuals,
    std::vector<float> &means,
    std::vector<float> &stds,
    std::vector<int> &modes)
{
    // first pass: mean, mode
    int gt_counts[3] = {0};
    size_t gt_sum = 0;
    size_t nan_count = 0;

    for (size_t i = 0; i < bedcol.size(); ++i)
    {
        unsigned char bed_byte = bedcol[i];
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); ++j)
        {
            size_t ix = (4 * bed_byte) + j;
            if (bed_lut_b[ix] == 0.0)
            {
                nan_count += 1;
            }
            else
            {
                gt_counts[(size_t)bed_lut_a[ix]] += 1;
                gt_sum += bed_lut_a[ix];
            }
        }
    }

    // find mode
    size_t mode = 0;
    for (size_t i = 1; i < 3; ++i)
    {
        if (gt_counts[i] > gt_counts[mode])
        {
            mode = i;
        }
    }

    float mean = gt_sum / (float)(num_individuals - nan_count);

    // second pass: impute, std
    float sum_squares = 0.0;
    for (size_t i = 0; i < bedcol.size(); ++i)
    {
        unsigned char bed_byte = bedcol[i];
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); ++j)
        {
            size_t ix = (4 * bed_byte) + j;
            if (bed_lut_b[ix] == 1.0)
            {
                sum_squares += (bed_lut_a[ix] - mean) * (bed_lut_a[ix] - mean);
            }
        }
    }

    float std = std::sqrt(sum_squares / (float)(num_individuals - nan_count));
    means.push_back(mean);
    stds.push_back(std);
    modes.push_back(mode);
}

void compute_bed_col_stats_impute(
    const std::vector<unsigned char> &bedcol,
    const size_t num_individuals,
    std::vector<float> &means,
    std::vector<float> &stds,
    std::vector<int> &modes)
{
    // first pass: mean, mode
    int gt_counts[3] = {0};
    size_t gt_sum = 0;
    size_t nan_count = 0;

    for (size_t i = 0; i < bedcol.size(); ++i)
    {
        unsigned char bed_byte = bedcol[i];
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); ++j)
        {
            size_t ix = (4 * bed_byte) + j;
            if (bed_lut_b[ix] == 0.0)
            {
                nan_count += 1;
            }
            else
            {
                gt_counts[(size_t)bed_lut_a[ix]] += 1;
                gt_sum += bed_lut_a[ix];
            }
        }
    }

    // find mode
    size_t mode = 0;
    for (size_t i = 1; i < 3; ++i)
    {
        if (gt_counts[i] > gt_counts[mode])
        {
            mode = i;
        }
    }

    gt_sum += (float)(mode * nan_count);
    float mean = gt_sum / (float)num_individuals;

    // second pass: impute, std
    float sum_squares = 0.0;
    for (size_t i = 0; i < bedcol.size(); ++i)
    {
        unsigned char new_byte = 0;
        unsigned char bed_byte = bedcol[i];
        for (size_t j = 0; j < 4; ++j)
        {
            new_byte <<= 2;
            size_t ix = (4 * bed_byte) + 3 - j;
            size_t curr_val = 2;
            if ((i * 4 + 3 - j) < num_individuals)
            {
                if (bed_lut_b[ix] == 0)
                {
                    curr_val = mode;
                }
                else
                {
                    curr_val = (size_t)bed_lut_a[ix];
                }
                float curr_diff = (float)curr_val - mean;
                sum_squares += curr_diff * curr_diff;
            }
            new_byte += gt_to_bed_value[curr_val];
        }
    }

    float std = std::sqrt(sum_squares / (float)num_individuals);
    means.push_back(mean);
    stds.push_back(std);
    modes.push_back(mode);
}

void prep_bed_no_impute(
    BfilesBase bfiles,
    std::string out_dir)
{
    if (!bfiles.has_valid_bed_prefix())
    {
        std::cout << "Invalid prefix bytes in bed" << std::endl;
        exit(1);
    }
    size_t num_individuals = count_lines(bfiles.fam());
    BimInfo bim(bfiles.bim());
    BedDims dim(num_individuals, bim.number_of_lines);
    size_t bed_block_size = dim.bytes_per_col();

    std::vector<unsigned char> bedcol(bed_block_size);

    std::vector<float> means;
    std::vector<float> stds;
    std::vector<int> modes;
    std::ifstream bed_file(bfiles.bed(), std::ios::binary);
    int gt_counts[3] = {0};

    while (bed_file.read(reinterpret_cast<char *>(bedcol.data()), bed_block_size))
    {
        compute_bed_col_stats_no_impute(bedcol, num_individuals, means, stds, modes);
    }

    dim.to_file(bfiles.dim());
    write_single_column_file_with_suffix(means, bfiles.means());
    write_single_column_file_with_suffix(stds, bfiles.stds());
    write_single_column_file_with_suffix(modes, bfiles.modes());
}