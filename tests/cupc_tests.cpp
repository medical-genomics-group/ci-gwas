#include <vector>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <mps/cuPC-S-sparse.h>
#include <mps/cu_test_helpers.h>
#include <mps/gpuerrors.h>
#include <test_data/bed_marker_test_set.h>

const int NUMBER_OF_LEVELS = 50;

auto std_normal_qnorm(const double p) -> double
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

auto threshold_array(const int n, const double alpha) -> std::array<double, NUMBER_OF_LEVELS>
{
    std::array<double, NUMBER_OF_LEVELS> thr{0.0};
    // my loop range is exclusive of the last i, so I don't subtract one here
    const int n_thr = (NUMBER_OF_LEVELS < (n - 3)) ? NUMBER_OF_LEVELS : n;
    const double half = 0.5;
    for (size_t i = 0; i < n_thr; i++)
    {
        thr[i] = abs(std_normal_qnorm(half * alpha) / sqrt(n - i - 3));
    }
    return thr;
}

TEST(cuPCSparseTests, CalIndepL1SingleBlock)
{
    const double alpha = 0.05;
    int w = 3;
    int p = BMT2_NUM_PHEN;
    int m = BMT2_NUM_MARKERS;
    int max_marker_degree = 2 * w + p;
    int max_phen_degree = m + p;
    size_t mixed_matrix_size = max_marker_degree * m + max_phen_degree * p;
    std::vector<int> G(mixed_matrix_size, 0);
    std::array<double, NUMBER_OF_LEVELS> Th = threshold_array(p + m, alpha);

    test_cal_Indepl0(
        bmt2_sparse_corrs,
        &m,
        &p,
        &w,
        G.data(),
        Th.data());
}
