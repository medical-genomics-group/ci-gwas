#include <gtest/gtest.h>
#include <mps/hetcor-cuPC-S.h>
#include <mps/cuPC_call_prep.h>
#include <test_data/cupc_test_set.h>

#include <cmath>
#include <vector>

TEST(hetcor_cuPC, expected_skeleton_n10)
{
    // call cuPC
    int max_level = 14;
    int p = N_N10;
    float threshold = hetcor_threshold(ALPHA_N10);
    std::vector<float> N(p * p, SAMPLE_SIZE_N10);
    const size_t sepset_size = p * p * ML;
    const size_t g_size = p * p;
    std::vector<int> G(g_size, 1);
    int l = 0;
    Skeleton(C_N10.data(), &p, G.data(), N.data(), &th, &l, &max_level);

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