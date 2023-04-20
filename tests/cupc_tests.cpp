#include <gtest/gtest.h>
#include <mps/cuPC-S.h>
#include <mps/cuPC_call_prep.h>
#include <test_data/cupc_test_set.h>

#include <cmath>
#include <vector>

TEST(threshold, at_10_e_min_8)
{
    // call cuPC
    int max_level = 14;
    std::vector<float> Th = threshold_array(500000, 0.00000001);
    for (size_t i = 0; i < max_level + 1; ++i)
    {
        EXPECT_EQ(Th[i], 0.007);
    }
}

TEST(cuPC, expected_skeleton_n10)
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
