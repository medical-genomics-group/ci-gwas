#include <gtest/gtest.h>
#include <mps/cuPC_call_prep.h>
#include <mps/cu_test_helpers.h>
#include <test_data/cupc_test_set.h>

#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <vector>

TEST(cuPCTests, ExpectedSkeletonN10)
{
    // call cuPC
    int max_level = 14;
    std::vector<float> Th = threshold_array(SAMPLE_SIZE_N10, ALPHA_N10);
    int p = N_N10;
    const size_t sepset_size = p * p * ML;
    const size_t g_size = p * p;
    std::vector<float> pmax(g_size, 0.0);
    std::vector<int> G(g_size, 1);
    std::vector<int> sepset(sepset_size, 0);
    int l = 0;
    Skeleton(C_N10.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data());

    // printf("ix | obs | exp \n");
    // for (size_t i = 0; i < CUPCT1_ADJSIZE; ++i)
    // {
    //     printf("%i | %i | %i \n", i, G[i], cupct1_g[i]);
    // }

    for (size_t i = 0; i < p * p; ++i)
    {
        EXPECT_EQ(G[i], A_N10[i]);
    }
}
