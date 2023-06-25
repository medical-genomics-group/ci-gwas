#include <mps/io.h>
#include <mps/trait_summary_stats.h>

TraitSummaryStats::TraitSummaryStats(const std::string path)
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
    std::cout << header << std::endl;
    num_phen = header.size();
    std::cout << num_phen << std::endl;
    corrs = std::vector(num_phen * num_phen, (float)0.0);
    int curr_corr_line = 0;

    while (std::getline(corr_file, line))
    {
        std::vector<std::string> fields = split_line(line);
        for (size_t j = 1; j <= num_phen; j++)
        {
            corrs[curr_corr_line * num_phen + j - 1] = std::stof(fields[j]);
        }
        curr_corr_line += 1;
    }

    // make corrs symmetric
    for (size_t i = 0; i < num_phen; i++)
    {
        for (size_t j = i + 1; j < num_phen; j++)
        {
            corrs[j * num_phen + i] = corrs[i * num_phen + j];
        }
    }
}
