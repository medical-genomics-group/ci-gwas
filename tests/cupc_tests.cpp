#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <mps/cu_test_helpers.h>
#include <test_data/cupc_test_set.h>
#include <vector>

const int NUMBER_OF_LEVELS = 50;

auto std_normal_qnorm(const float p) -> float
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

auto threshold_array(const int n, const float alpha) -> std::array<float, NUMBER_OF_LEVELS>
{
    std::array<float, NUMBER_OF_LEVELS> thr{0.0};
    // my loop range is exclusive of the last i, so I don't subtract one here
    const int n_thr = (NUMBER_OF_LEVELS < (n - 3)) ? NUMBER_OF_LEVELS : n;
    const float half = 0.5;
    for (size_t i = 0; i < n_thr; i++)
    {
        thr[i] = abs(std_normal_qnorm(half * alpha) / sqrt(n - i - 3));
    }
    return thr;
}

TEST(cuPCSparseTests, CalIndepL0SingleBlock)
{
    const float alpha = 0.1;
    int w = CUPCT1_W;
    int p = CUPCT1_P;
    int m = CUPCT1_M;
    int max_marker_degree = 2 * w + p;
    int max_phen_degree = m + p;
    size_t mixed_matrix_size = max_marker_degree * m + max_phen_degree * p;
    std::vector<int> G(mixed_matrix_size, 0);
    std::array<float, NUMBER_OF_LEVELS> Th = threshold_array(m + p, alpha);

    test_cal_Indepl0(
        &cupct1_c[0],
        &m,
        &p,
        &w,
        G.data(),
        Th.data());

    // printf("ix | obs | exp \n");
    // for (size_t i = 0; i < CUPCT1_ADJSIZE; ++i)
    // {
    //     printf("%i | %i | %i \n", i, G[i], cupct1_g[i]);
    // }

    for (size_t i = 0; i < CUPCT1_ADJSIZE; ++i)
    {
        EXPECT_EQ(G[i], cupct1_g[i]);
    }
}

TEST(cuPCSparseTests, CalIndepL0MultiBlock)
{
    const float alpha = 0.1;
    int w = CUPCT2_W;
    int p = CUPCT2_P;
    int m = CUPCT2_M;
    int max_marker_degree = 2 * w + p;
    int max_phen_degree = m + p;
    size_t mixed_matrix_size = max_marker_degree * m + max_phen_degree * p;
    std::vector<int> G(mixed_matrix_size, 0);
    std::array<float, NUMBER_OF_LEVELS> Th = threshold_array(m + p, alpha);

    test_cal_Indepl0(
        &cupct2_c[0],
        &m,
        &p,
        &w,
        G.data(),
        Th.data());

    printf("ix | obs | exp \n");
    for (size_t i = 0; i < CUPCT2_ADJSIZE; ++i)
    {
        if (G[i] != cupct2_g[i]) {
            printf("%i | %i | %i \n", i, G[i], cupct2_g[i]);
        }
    }

    for (size_t i = 0; i < CUPCT2_ADJSIZE; ++i)
    {
        EXPECT_EQ(G[i], cupct2_g[i]);
    }
}
