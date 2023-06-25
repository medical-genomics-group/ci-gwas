#include <mps/io.h>
#include <mps/marker_summary_stats.h>

MarkerSummaryStats::MarkerSummaryStats(const std::string path)
{
    std::vector<float> triu_corrs = read_floats_from_binary(path);
    num_markers = ((sqrt(4 * 2 * triu_corrs.size() + 1) - 1) / 2) + 1;
    corrs = std::vector(num_markers * num_markers, (float)0.0);
    size_t corr_ix = 0;
    for (size_t i = 0; i < num_markers; i++)
    {
        for (size_t j = i + 1; j < num_markers; j++)
        {
            corrs[i * num_markers + j] = triu_corrs[corr_ix];
            corrs[j * num_markers + i] = triu_corrs[corr_ix];
            corr_ix += 1;
        }
    }
}
