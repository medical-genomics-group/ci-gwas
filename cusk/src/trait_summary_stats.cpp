#include <mps/io.h>
#include <math.h>
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
    num_phen = header.size();
    corrs = std::vector(num_phen * num_phen, (float)1.0);
    int curr_corr_line = 0;

    while (std::getline(corr_file, line))
    {
        std::vector<std::string> fields = split_line(line);

        if (fields.empty())
        {
            break;
        }

        for (size_t j = 1; j <= num_phen; j++)
        {
            float corr_val = zero_if_nan(std::stof(fields[j]));
            corrs[curr_corr_line * num_phen + j - 1] = corr_val;
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

TraitSummaryStats::TraitSummaryStats(
    const std::string corr_path,
    const std::string se_path)
{
    float corr_val;
    float sample_size;
    float se;
    float ss_sqrt;
    std::string corr_line;
    std::string se_line;
    std::ifstream corr_file(corr_path);
    std::ifstream se_file(se_path);

    // read header
    if (!(std::getline(corr_file, corr_line)))
    {
        std::cout << "trait summary corr file seems to be empty" << std::endl;
        exit(1);
    }
    // read header
    if (!(std::getline(se_file, se_line)))
    {
        std::cout << "trait summary se file seems to be empty" << std::endl;
        exit(1);
    }

    /*
    I am assuming here that the corr file has the same shape as the se file,
    without implementing any checks:
    TODO: implement checks, e.g. same number of fields, same number of lines
    */

    header = split_line(corr_line);
    num_phen = header.size();
    corrs = std::vector(num_phen * num_phen, (float)1.0);
    sample_sizes = std::vector(num_phen * num_phen, (float)0.0);
    int curr_corr_line = 0;

    while (std::getline(corr_file, corr_line))
    {
        std::getline(se_file, se_line);
        std::vector<std::string> fields_corr = split_line(corr_line);
        std::vector<std::string> fields_se = split_line(se_line);

        if (fields_corr.empty())
        {
            break;
        }

        for (size_t j = 1; j <= num_phen; j++)
        {
            corr_val = std::stof(fields_corr[j]);
            if (std::isnan(corr_val)) {
                corr_val = 0.0;
                sample_size = NAN;
            } else {
                se = std::stof(fields_se[j]);
                ss_sqrt = (1.0 - (corr_val * corr_val)) / se;
                sample_size = ss_sqrt * ss_sqrt;
            }
            corrs[curr_corr_line * num_phen + j - 1] = corr_val;
            sample_sizes[curr_corr_line * num_phen + j - 1] = sample_size;
        }
        curr_corr_line += 1;
    }

    // make matrices symmetric
    for (size_t i = 0; i < num_phen; i++)
    {
        for (size_t j = i + 1; j < num_phen; j++)
        {
            corrs[j * num_phen + i] = corrs[i * num_phen + j];
            sample_sizes[j * num_phen + i] = sample_sizes[i * num_phen + j];
        }
    }
}
