#include <mps/io.h>
#include <mps/marker_summary_stats.h>

/**
 * @brief Load mxm matrix from binary file. Assumed to contain float32 values in lower triangular
 * matrix format, including the diagonal, row major.
 */
MarkerSummaryStats::MarkerSummaryStats(const std::string path)
{
    std::vector<float> tril_corrs = read_floats_from_binary(path);
    num_markers = ((sqrt(4 * 2 * tril_corrs.size() + 1) - 1) / 2);
    corrs = std::vector(num_markers * num_markers, (float)1.0);
    size_t corr_ix = 0;
    for (size_t i = 0; i < num_markers; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            float corr_val = zero_if_nan(tril_corrs[corr_ix]);
            corrs[i * num_markers + j] = corr_val;
            corrs[j * num_markers + i] = corr_val;
            corr_ix += 1;
        }
    }
}
