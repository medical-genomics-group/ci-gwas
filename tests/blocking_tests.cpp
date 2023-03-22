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

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(block_chr, expected_results_chr1)
{
    std::vector<MarkerBlock> obs = block_chr(CHR_ONE_MCORRK_SUM_ABS_ROWS, "1", 11000);
    std::vector<MarkerBlock> exp = {
        MarkerBlock("1", 0, 5090),        MarkerBlock("1", 5091, 11133),
        MarkerBlock("1", 11134, 15892),   MarkerBlock("1", 15893, 22285),
        MarkerBlock("1", 22286, 27105),   MarkerBlock("1", 27106, 30102),
        MarkerBlock("1", 30103, 37499),   MarkerBlock("1", 37500, 42150),
        MarkerBlock("1", 42151, 45358),   MarkerBlock("1", 45359, 49723),
        MarkerBlock("1", 49724, 58066),   MarkerBlock("1", 58067, 63058),
        MarkerBlock("1", 63059, 67599),   MarkerBlock("1", 67600, 75794),
        MarkerBlock("1", 75795, 81929),   MarkerBlock("1", 81930, 87765),
        MarkerBlock("1", 87766, 92475),   MarkerBlock("1", 92476, 95480),
        MarkerBlock("1", 95481, 99014),   MarkerBlock("1", 99015, 102591),
        MarkerBlock("1", 102592, 110967), MarkerBlock("1", 110968, 116784),
        MarkerBlock("1", 116785, 123611), MarkerBlock("1", 123612, 129527),
        MarkerBlock("1", 129528, 139781), MarkerBlock("1", 139782, 144870),
        MarkerBlock("1", 144871, 151294), MarkerBlock("1", 151295, 157216),
        MarkerBlock("1", 157217, 162050), MarkerBlock("1", 162051, 165281),
        MarkerBlock("1", 165282, 169903)};

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}
