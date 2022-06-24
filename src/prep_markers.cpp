// input: .bed path, .bim path, .fam path
// 
// we want to:
// - split .bed by chromosome
// - standardize columns of X (blocks in the .bed file)
// - compute means and standard deviations of columns (blocks in the .bed file)
//
// we need number of individuals and number of markers.
// the number of markers P is equal to the number of rows in the .bim file.
// the number of individuals N is equal to the number of rows in .fam file.

#include <cmath>
#include <iterator>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <mps/prep_markers.h>
#include <mps/bed_lut.h>

auto count_lines(std::string file_path) -> int {
    int number_of_lines = 0;
    std::string line;
    std::ifstream fin(file_path);

    while (std::getline(fin, line))
        ++number_of_lines;

    return number_of_lines;
}

void split_bim_line(std::string line, std::string *buf)
{
    std::istringstream ss(line);
 
    std::string word;
    size_t word_ix = 0;
    
    while (ss >> word)
    {
        assert((word_ix < 6) && "found more than 6 columns in bim file, aborting.");
        buf[word_ix] = word;
        word_ix += 1;
    }
}

auto parse_bim(std::string bim_path) -> bimInfo {
    bimInfo res;
    res.number_of_lines = 0;
    res.chr_ids = {};
    res.num_markers_on_chr = {};

    std::string bim_line[6];
    std::string line;
    std::ifstream bim(bim_path);

    while (std::getline(bim, line)) {
        split_bim_line(line, bim_line);
        if ((res.number_of_lines == 0) || (bim_line[0] != res.chr_ids.back())) {
            res.chr_ids.push_back(bim_line[0]);
            res.num_markers_on_chr.push_back(0);
        }
        ++res.num_markers_on_chr.back();
        ++res.number_of_lines;
    }

    return res;
}

auto make_path(std::string out_dir, std::string chr_id, std::string suffix) -> std::string {
    std::string filename = chr_id + suffix;
    std::string outpath;
    if (out_dir.back() != '/') {
        outpath = out_dir + '/' + filename;
    } else {
        outpath = out_dir + filename;
    }
    return outpath;
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
        fout << chr_marker_means[i];
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
        fout << chr_marker_stds[i];
    }

    fout.close();       
}


// Split .bed by chromosome, impute NaN to median, compute col means and std
void prep_bed(std::string bed_path,
              std::string bim_path,
              std::string fam_path,
              std::string out_dir,
              size_t mem_gb) // max memory used by this fn)
{
    size_t num_individuals = count_lines(fam_path);
    bimInfo bim_info =  parse_bim(bim_path);
    size_t bed_block_size = (num_individuals + 3) / 4;

    std::vector<unsigned char> bed_block(bed_block_size);
    std::ifstream bed_file(bed_path, std::ios::binary);

    // read first three bytes, check if correct code
    bed_file.read(reinterpret_cast<char*>(bed_block.data()), 3);
    unsigned char bed_2_code[3] = {0x6c, 0x1b, 0x00};
    for (size_t i = 0; i < 3; ++i) {
        assert((bed_block[i] == bed_2_code[i]) && "unexpected magic number in bed file.");
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

    while (bed_file.read(reinterpret_cast<char*>(bed_block.data()), bed_block_size)) {
        // first pass: mean, median
        gt_counts[0] = 0;
        gt_counts[1] = 0;
        gt_counts[2] = 0;
        size_t gt_sum = 0;
        size_t nan_count = 0;

        for (size_t i = 0; i < bed_block_size; ++i) {
            unsigned char bed_byte = bed_block[i];
            for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); ++j) {
                size_t ix = bed_byte + j;
                if (bed_lut_b[ix] == 0.0) {
                    nan_count += 1;
                } else {
                    gt_counts[(size_t)bed_lut_a[ix]] += 1;
                    gt_sum += bed_lut_a[ix];
                }
            }
        }
        
        // find median
        size_t median = 0;
        for (size_t i = 1; i < 3; ++i) {
            if (gt_counts[i] > gt_counts[median]) {
                median = i;
            }
        }

        gt_sum += (float)(median * nan_count);
        float mean = gt_sum / (float)num_individuals;

        // second pass: impute, std
        float sum_squares = 0.0;
        for (size_t i = 0; i < bed_block_size; ++i) {
            unsigned char new_byte = 0;
            unsigned char bed_byte = bed_block[i];
            for (size_t j = 0; j < 4; ++j) {
                size_t ix = bed_byte + 3 - j;
                if ((i * 4 + j) < num_individuals) {
                    new_byte <<= 2;
                    size_t curr_val;
                    if (bed_lut_b[ix] == 0) {
                        curr_val = median;
                    } else {
                        curr_val = (size_t)bed_lut_a[ix];
                    }
                    new_byte += gt_to_bed_value[curr_val];
                    float curr_diff = (float)curr_val - mean;
                    sum_squares += curr_diff * curr_diff;
                }
            }
            out_buf.push_back(new_byte);
        }

        float std = std::sqrt(sum_squares / (float)num_individuals);
        chr_marker_means.push_back(mean);
        chr_marker_stds.push_back(std);
        chr_processed_markers += 1;

        // write to file
        if (
                (blocks_in_buf >= max_num_blocks) || 
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
            if (chr_processed_markers == bim_info.num_markers_on_chr[chr_ix]) {
                chr_ix += 1;
                chr_processed_markers = 0;
            }
        }
    }
}


