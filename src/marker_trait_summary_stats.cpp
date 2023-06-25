#include <mps/io.h>
#include <mps/marker_trait_summary_stats.h>

#include <string>

MarkerTraitSummaryStats::MarkerTraitSummaryStats(const std::string path)
{
    std::string line;
    std::ifstream corr_file(path);

    // read header
    if (!(std::getline(corr_file, line)))
    {
        std::cout << "marker-trait summary stat file seems to be empty" << std::endl;
        exit(1);
    }

    header = split_line(line);

    if ((header[0] != "chr") || (header[1] != "snp") || (header[2] != "ref"))
    {
        std::cout << "marker-trait summary stat file has bad header" << std::endl;
        exit(1);
    }

    num_phen = header.size() - 3;
    num_markers = 0;
    corrs = {};

    while (std::getline(corr_file, line))
    {
        split_line(line);
        for (size_t j = 3; j < num_phen + 3; j++)
        {
            corrs.push_back(std::stof(line[j]));
        }
        num_markers += 1;
    }
}
