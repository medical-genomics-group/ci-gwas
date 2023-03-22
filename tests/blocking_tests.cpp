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

/*
TEST(block_chr, expected_results_chr1)
{
    std::vector<MarkerBlock> obs = block_chr(CHR_ONE_MCORRK_SUM_ABS_ROWS, "1", 11000);
    std::vector<MarkerBlock> exp = {
        MarkerBlock("1", 0, 5087),        MarkerBlock("1", 5088, 11133),
        MarkerBlock("1", 11134, 15892),   MarkerBlock("1", 15893, 22285),
        MarkerBlock("1", 22286, 27105),   MarkerBlock("1", 27106, 30102),
        MarkerBlock("1", 30103, 37499),   MarkerBlock("1", 37500, 42150),
        MarkerBlock("1", 42151, 45358),   MarkerBlock("1", 45359, 49722),
        MarkerBlock("1", 49723, 58066),   MarkerBlock("1", 58067, 63057),
        MarkerBlock("1", 63058, 67598),   MarkerBlock("1", 67599, 75794),
        MarkerBlock("1", 75795, 81929),   MarkerBlock("1", 81930, 87764),
        MarkerBlock("1", 87765, 92475),   MarkerBlock("1", 92476, 95480),
        MarkerBlock("1", 95481, 99015),   MarkerBlock("1", 99016, 102590),
        MarkerBlock("1", 102591, 110967), MarkerBlock("1", 110968, 116784),
        MarkerBlock("1", 116785, 123611), MarkerBlock("1", 123612, 129526),
        MarkerBlock("1", 129527, 139781), MarkerBlock("1", 139782, 144869),
        MarkerBlock("1", 144870, 151293), MarkerBlock("1", 151294, 157216),
        MarkerBlock("1", 157217, 162049), MarkerBlock("1", 162050, 165281),
        MarkerBlock("1", 165282, 169903)};

    ASSERT_EQ(exp.size(), obs.size());

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}
*/

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
