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
        MarkerBlock("1", 0, 194),
        MarkerBlock("1", 195, 335),
        MarkerBlock("1", 336, 620),
        MarkerBlock("1", 621, 843),
        MarkerBlock("1", 844, 1227),
        MarkerBlock("1", 1228, 1447),
        MarkerBlock("1", 1448, 1910),
        MarkerBlock("1", 1911, 2112),
        MarkerBlock("1", 2113, 2504),
        MarkerBlock("1", 2505, 2735),
        MarkerBlock("1", 2736, 2930),
        MarkerBlock("1", 2931, 3085),
        MarkerBlock("1", 3086, 3172),
        MarkerBlock("1", 3173, 3352),
        MarkerBlock("1", 3353, 3574),
        MarkerBlock("1", 3575, 3897),
        MarkerBlock("1", 3898, 3997)};

    EXPECT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(block_chr, expected_results_chr1)
{
    std::vector<MarkerBlock> obs = block_chr(CHR_ONE_MCORRK_SUM_ABS_ROWS, "1", 200);
    std::vector<MarkerBlock> exp = {
        MarkerBlock("1", 0, 84),
        MarkerBlock("1", 85, 151),
        MarkerBlock("1", 152, 181),
        MarkerBlock("1", 182, 283),
        MarkerBlock("1", 284, 351),
        MarkerBlock("1", 352, 403),
        MarkerBlock("1", 404, 469),
        MarkerBlock("1", 470, 563),
        MarkerBlock("1", 564, 710),
        MarkerBlock("1", 711, 787),
        MarkerBlock("1", 788, 819),
        MarkerBlock("1", 820, 903),
        MarkerBlock("1", 904, 999)};

    ASSERT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(hanning_smoothing, expected_results)
{
    std::vector<float> obs =
        hanning_smoothing(std::vector<float>(TEST_V.begin(), TEST_V.begin() + 1000), 101);
    std::vector<float> exp = TEST_V_SMOOTH;

    ASSERT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_NEAR(exp[i], obs[i], 0.01);
    }
}
