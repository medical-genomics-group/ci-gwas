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

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <mps/prep_markers.h>

auto count_lines(std::string file_path) -> int {
    int number_of_lines = 0;
    std::string line;
    std::ifstream fin(file_path);

    while (std::getline(fin, line))
        ++number_of_lines;

    return number_of_lines;
}

void split_line(std::string line, std::string *buf)
{
    std::istringstream ss(line);
 
    std::string word;
    size_t word_ix = 0;
    
    while (ss >> word)
    {
        buf[word_ix] = word;
        word_ix += 1;
    }
}

auto parse_bim(std::string bim_path) -> bimInfo {
    bimInfo res;
    res.number_of_lines = 0;
    res.chr_ids = {};
    res.chr_first_ix = {};

    std::string bim_line[6];
    std::string line;
    std::ifstream bim(bim_path);

    while (std::getline(bim, line)) {
        split_line(line, bim_line);
        if ((res.number_of_lines == 0) || (bim_line[0] != res.chr_ids.back())) {
            res.chr_ids.push_back(bim_line[0]);
            res.chr_first_ix.push_back(res.number_of_lines);
        }
        ++res.number_of_lines;
    }

    return res;
}
