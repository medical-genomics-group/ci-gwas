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
        if (line_num > block.get_last_marker_global_ix())
        {
            break;
        }
        else if (line_num >= block.get_first_marker_global_ix())
        {
            std::vector<std::string> fields = split_line(line);
            for (size_t j = 3; j < num_phen + 3; j++)
            {
                if ((fields[j] == "NA") || (fields[j] == "NaN") || (fields[j] == "nan"))
                {
                    corrs.push_back(0.0);
                }
                else
                {
                    corrs.push_back(std::stof(fields[j]));
                }
            }
            num_markers += 1;
        }
        line_num += 1;
    }
}

MarkerTraitSummaryStats::MarkerTraitSummaryStats(
    const std::string corr_path,
    const std::string se_path,
    const MarkerBlock block,
    const float nan_sample_size_fill_val)
{
    std::string line;
    std::ifstream corr_file(corr_path);

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
        if (line_num > block.get_last_marker_global_ix())
        {
            break;
        }
        else if (line_num >= block.get_first_marker_global_ix())
        {
            std::vector<std::string> fields = split_line(line);
            for (size_t j = 3; j < num_phen + 3; j++)
            {
                if ((fields[j] == "NA") || (fields[j] == "NaN") || (fields[j] == "nan"))
                {
                    corrs.push_back(0.0);
                }
                else
                {
                    corrs.push_back(std::stof(fields[j]));
                }
            }
            num_markers += 1;
        }
        line_num += 1;
    }

    std::string line;
    std::ifstream se_file(se_path);

    // read header
    if (!(std::getline(se_file, line)))
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
    sample_sizes = {};

    size_t line_num = 0;
    while (std::getline(se_file, line))
    {
        if (line_num > block.get_last_marker_global_ix())
        {
            break;
        }
        else if (line_num >= block.get_first_marker_global_ix())
        {
            std::vector<std::string> fields = split_line(line);
            for (size_t j = 3; j < num_phen + 3; j++)
            {
                if ((fields[j] == "NA") || (fields[j] == "NaN") || (fields[j] == "nan"))
                {
                    sample_sizes.push_back(nan_sample_size_fill_val);
                }
                else
                {
                    sample_sizes.push_back(std::stof(fields[j]));
                }
            }
            num_markers += 1;
        }
        line_num += 1;
    }
}

MarkerTraitSummaryStats::MarkerTraitSummaryStats(
    const std::string path, const std::vector<int> marker_ixs
)
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
        if (line_num == marker_ixs[num_markers])
        {
            std::vector<std::string> fields = split_line(line);
            for (size_t j = 3; j < num_phen + 3; j++)
            {
                if ((fields[j] == "NA") || (fields[j] == "NaN") || (fields[j] == "nan"))
                {
                    corrs.push_back(0.0);
                }
                else
                {
                    corrs.push_back(std::stof(fields[j]));
                }
            }
            num_markers += 1;
        }
        line_num += 1;
    }
}