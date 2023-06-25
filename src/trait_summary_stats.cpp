#include <mps/io.h>
#include <mps/trait_summary_stats.h>

TraitSummaryStats(const std::string path)
{
    std::string line;
    std::ifstream corr_file(path);

    // read header
    if (!(std::getline(corr_file, line)))
    {
        std::cout << "trait summary stat file seems to be empty" << std::endl;
        exit(1);
    }

    header = split_line(line);
    num_phen = header.size();
    corrs = {};

    while (std::getline(corr_file, line))
    {
        split_line(line);
        for (size_t i = 1; i <= num_phen; i++)
        {
            corrs.push_back(std::stof(line[i]));
        }
    }

    // make corrs symmetric
    for (size_t i = 0; i < num_phen; i++)
    {
        for (size_t j = i + 1; j < num_phen; j++)
        {
            corrs[j * num_phen + i] = corrs[i * num_phen + j]
        }
    }
}
