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
        MarkerBlock("1", 0, 249),
        MarkerBlock("1", 250, 537),
        MarkerBlock("1", 538, 803),
        MarkerBlock("1", 804, 1138),
        MarkerBlock("1", 1139, 1406),
        MarkerBlock("1", 1407, 1553),
        MarkerBlock("1", 1554, 1830),
        MarkerBlock("1", 1831, 2100),
        MarkerBlock("1", 2101, 2506),
        MarkerBlock("1", 2507, 2723),
        MarkerBlock("1", 2724, 2934),
        MarkerBlock("1", 2935, 3145),
        MarkerBlock("1", 3146, 3358),
        MarkerBlock("1", 3359, 3571),
        MarkerBlock("1", 3572, 3898),
        MarkerBlock("1", 3899, 3998)};

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}
