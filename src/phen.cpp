#include <fstream>
#include <iostream>
#include <mps/io.h>
#include <mps/phen.h>

// Return only the phenotype value part of each line
auto split_phen_line(std::string line) -> std::vector<float>
{
    std::vector<std::string> full_line = split_line(line);
    std::vector<float> res;
    for (size_t i = 2; i < full_line.size(); ++i)
    {
        res.push_back((float)atof(full_line[i].c_str()));
    }
    res;
}

auto load_phen(std::string phen_path) -> Phen
{
    Phen res;
    std::string line;
    std::ifstream phen(phen_path);
    std::vector<float> row_major_vals;

    while (std::getline(phen, line))
    {
        std::vector<float> phen_vals = split_phen_line(line);
        if (res.num_phen == 0)
        {
            res.num_phen = phen_vals.size();
        }
        else if (res.num_phen != phen_vals.size())
        {
            printf("Inconsistent row width in .phen file");
            exit(1);
        }
        for (float v : phen_vals)
        {
            row_major_vals.push_back(v);
        }
        ++res.num_samples;
    }
    // make col major
    for (size_t p = 0; p < res.num_phen; ++p)
    {
        for (size_t i = 0; i < res.num_samples; ++i)
        {
            res.data.push_back(row_major_vals[(res.num_phen * i) + p]);
        }
    }
    return res;
}