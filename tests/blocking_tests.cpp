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
        MarkerBlock("1", 0, 147),
        MarkerBlock("1", 148, 249),
        MarkerBlock("1", 250, 585),
        MarkerBlock("1", 586, 763),
        MarkerBlock("1", 764, 1259),
        MarkerBlock("1", 1260, 1570),
        MarkerBlock("1", 1571, 1680),
        MarkerBlock("1", 1681, 1785),
        MarkerBlock("1", 1786, 1999)};

    for (size_t i = 0; i < exp.size(); ++i)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}
