#include <math.h>
#include <mps/io.h>

#include <vector>

const int MAX_BLOCK_SIZE_TOL = 100;

double hann(const int n, const int m)
{
    return 0.5 - 0.5 * cosf(2.0 * M_PI * (double)n / ((double)m - 1.0));
}

std::vector<double> hanning_smoothing(const std::vector<float> &v, const int window_size)
{
    std::vector<double> hanning_window;
    for (int i = 0; i < window_size; i++)
    {
        hanning_window.push_back(hann(i, window_size));
    }

    int margin = window_size / 2;

    std::vector<double> res(v.size(), 0.0);

    for (int center = margin; center < (v.size() - margin); center++)
    {
        for (int i = 0; i < window_size; i++)
        {
            res[center] += hanning_window[i] * (double)v[center - margin + i];
        }
    }

    return res;
}

std::vector<size_t> local_minima(const std::vector<double> &v)
{
    std::vector<size_t> res;
    double left = 0.0;
    for (size_t i = 1; i < (v.size() - 1); i++)
    {
        if ((left > v[i]) && (v[i] < v[i + 1]))
        {
            res.push_back(i);
            left = 0.0;
        }
        else if (v[i] > left)
        {
            left = v[i];
        }
    }
    return res;
}

std::vector<MarkerBlock> blocks_from_minima(
    const std::vector<size_t> &minima, const std::string chr_id, const size_t num_vars
)
{
    std::vector<MarkerBlock> res;
    size_t prev = 0;
    for (auto pos : minima)
    {
        res.push_back(MarkerBlock(chr_id, prev, pos, 0));
        prev = pos + 1;
    }
    res.push_back(MarkerBlock(chr_id, prev, num_vars - 1, 0));
    return res;
}

int largest_block_size(const std::vector<MarkerBlock> &blocks)
{
    size_t res = 0;
    for (auto block : blocks)
    {
        if (block.block_size() > res)
        {
            res = block.block_size();
        }
    }
    return res;
}

std::vector<MarkerBlock> block_chr_with_window_size(
    const std::vector<float> &forward_corr_sums, const std::string chr_id, const int window_size
)
{
    int n = forward_corr_sums.size();
    std::vector<double> smooth = hanning_smoothing(forward_corr_sums, window_size);
    std::vector<size_t> minima = local_minima(smooth);
    return blocks_from_minima(minima, chr_id, n);
}

int make_odd(int v)
{
    if (v % 2 == 0)
    {
        return v - 1;
    }
    return v;
}

std::vector<MarkerBlock> block_chr(
    const std::vector<float> &forward_corr_sums, const std::string chr_id, const int max_block_size
)
{
    int too_large = forward_corr_sums.size();
    int too_small = 3;
    int window_size = make_odd((too_large + too_small) / 2);

    std::vector<MarkerBlock> res =
        block_chr_with_window_size(forward_corr_sums, chr_id, window_size);
    int lbs = largest_block_size(res);

    while ((abs(lbs - max_block_size) > MAX_BLOCK_SIZE_TOL) || (lbs > max_block_size))
    {
        if (lbs > max_block_size)
        {
            too_large = std::min(too_large, window_size);
        }
        else
        {
            too_small = std::max(too_small, window_size);
        }
        int new_window_size = make_odd((too_large + too_small) / 2);
        if (new_window_size == window_size)
        {
            break;
        }
        window_size = new_window_size;

        res = block_chr_with_window_size(forward_corr_sums, chr_id, window_size);
        lbs = largest_block_size(res);
    }

    return res;
}
