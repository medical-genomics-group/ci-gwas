#include <mps/io.h>
#include <mps/marker_trait_summary_stats.h>

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
        std::vector<std::string> fields = split_line(line);
        for (size_t j = 3; j < num_phen + 3; j++)
        {
            corrs.push_back(std::stof(fields[j]));
        }
        num_markers += 1;
    }
}

MarkerTraitSummaryStats::MarkerTraitSummaryStats(const std::string path, const MarkerBlock block)
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

    size_t line_num = 0;
    while (std::getline(corr_file, line))
    {
        if (line_num > block.get_last_marker_ix())
        {
            break;
        }
        else if (line_num >= block.get_first_marker_ix())
        {
            std::vector<std::string> fields = split_line(line);
            for (size_t j = 3; j < num_phen + 3; j++)
            {
                std::cout << fields[j] << std::endl;
                corrs.push_back(std::stof(fields[j]));
            }
            num_markers += 1;
        }
        line_num += 1;
    }
}