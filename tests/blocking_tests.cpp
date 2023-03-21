#include <gtest/gtest.h>
#include <mps/blocking.h>
#include <test_data/blocking_test_set.h>

#include <filesystem>
#include <string>
#include <vector>

TEST(BlockChr, ExpectedResults)
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
