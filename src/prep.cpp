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

// Split .bed by chromosome, impute NaN to mode, compute col means and std
void prep_bed(
    std::string bed_path,
    std::string bim_path,
    std::string fam_path,
    std::string out_dir,
    size_t mem_gb) // max memory used by this fn)
{
    size_t num_individuals = count_lines(fam_path);
    bimInfo bim_info = prep_bim(bim_path, out_dir);
    size_t bed_block_size = (num_individuals + 3) / 4;

    for (size_t i = 0; i < bim_info.chr_ids.size(); ++i)
    {
        write_dims(num_individuals, bim_info.num_markers_on_chr[i], out_dir, bim_info.chr_ids[i]);
    }

    std::vector<unsigned char> bed_block(bed_block_size);
    std::ifstream bed_file(bed_path, std::ios::binary);

    // read first three bytes, check if correct code
    bed_file.read(reinterpret_cast<char *>(bed_block.data()), 3);
    unsigned char bed_col_maj_code[3] = {0x6c, 0x1b, 0x01};
    for (size_t i = 0; i < 3; ++i)
    {
        assert((bed_block[i] == bed_col_maj_code[i]) && "unexpected magic number in bed file.");
    }

    size_t max_num_bytes = mem_gb * 1'000'000'000;
    size_t max_num_blocks = (max_num_bytes / (bed_block_size + 2 * sizeof(float))) - 1;

    std::vector<float> chr_marker_means = {};
    std::vector<float> chr_marker_stds = {};
    std::vector<unsigned char> out_buf = {};
    out_buf.reserve(max_num_bytes);
    size_t blocks_in_buf = 0;
    size_t gt_counts[3] = {0};
    size_t chr_ix = 0;
    size_t chr_processed_markers = 0;

    while (bed_file.read(reinterpret_cast<char *>(bed_block.data()), bed_block_size))
    {
        // first pass: mean, mode
        gt_counts[0] = 0;
        gt_counts[1] = 0;
        gt_counts[2] = 0;
        size_t gt_sum = 0;
        size_t nan_count = 0;

        for (size_t i = 0; i < bed_block_size; ++i)
        {
            unsigned char bed_byte = bed_block[i];
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
        for (size_t i = 0; i < bed_block_size; ++i)
        {
            unsigned char new_byte = 0;
            unsigned char bed_byte = bed_block[i];
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
            out_buf.push_back(new_byte);
        }

        float std = std::sqrt(sum_squares / (float)num_individuals);
        chr_marker_means.push_back(mean);
        chr_marker_stds.push_back(std);
        chr_processed_markers += 1;

        // write to file
        if ((blocks_in_buf >= max_num_blocks) ||
            (chr_processed_markers == bim_info.num_markers_on_chr[chr_ix]))
        {
            write_bed(out_buf, out_dir, bim_info.chr_ids[chr_ix]);
            write_means(chr_marker_means, out_dir, bim_info.chr_ids[chr_ix]);
            write_stds(chr_marker_stds, out_dir, bim_info.chr_ids[chr_ix]);
            chr_marker_means.clear();
            chr_marker_stds.clear();
            // TODO: add magic numbers (maybe)
            out_buf.clear();
            blocks_in_buf = 0;
            if (chr_processed_markers == bim_info.num_markers_on_chr[chr_ix])
            {
                chr_ix += 1;
                chr_processed_markers = 0;
            }
        }
    }
}