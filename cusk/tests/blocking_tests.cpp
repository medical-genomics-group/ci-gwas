#include <gtest/gtest.h>
#include <mps/blocking.h>
#include <test_data/blocking_test_set.h>

#include <filesystem>
#include <string>
#include <vector>

TEST(block_chr, expected_results_synthetic_data)
{
    std::vector<MarkerBlock> obs = block_chr(TEST_V, "1", 500);
    std::vector<MarkerBlock> exp = {
        MarkerBlock("1", 0, 194, 0),
        MarkerBlock("1", 195, 335, 0),
        MarkerBlock("1", 336, 620, 0),
        MarkerBlock("1", 621, 843, 0),
        MarkerBlock("1", 844, 1227, 0),
        MarkerBlock("1", 1228, 1447, 0),
        MarkerBlock("1", 1448, 1910, 0),
        MarkerBlock("1", 1911, 2112, 0),
        MarkerBlock("1", 2113, 2504, 0),
        MarkerBlock("1", 2505, 2735, 0),
        MarkerBlock("1", 2736, 2930, 0),
        MarkerBlock("1", 2931, 3085, 0),
        MarkerBlock("1", 3086, 3172, 0),
        MarkerBlock("1", 3173, 3352, 0),
        MarkerBlock("1", 3353, 3574, 0),
        MarkerBlock("1", 3575, 3897, 0),
        MarkerBlock("1", 3898, 3997, 0),
    };

    EXPECT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(hanning_smoothing, expected_results)
{
    std::vector<double> obs =
        hanning_smoothing(std::vector<float>(TEST_V.begin(), TEST_V.begin() + 1000), 101);
    std::vector<double> exp = TEST_V_SMOOTH;

    ASSERT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_NEAR(exp[i], obs[i], 0.01);
    }
}
