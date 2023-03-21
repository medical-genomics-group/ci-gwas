#include <gtest/gtest.h>
#include <mps/parent_set.h>
#include <test_data/parent_set_test_set.h>

#include <vector>

TEST(parent_set, ExpectedResultsDepth1)
{
    std::vector<int> obs =
        parent_set(TEST_ADJ_MAT, TEST_NUM_MARKERS + TEST_NUM_PHEN, TEST_NUM_MARKERS, 1);
    for (size_t i = 0; i < obs.size(); ++i)
    {
        EXPECT_EQ(TEST_PAR_SET_D1[i], obs[i]);
    }
}

TEST(parent_set, ExpectedResultsDepth2)
{
    std::vector<int> obs =
        parent_set(TEST_ADJ_MAT, TEST_NUM_MARKERS + TEST_NUM_PHEN, TEST_NUM_MARKERS, 2);
    for (size_t i = 0; i < obs.size(); ++i)
    {
        EXPECT_EQ(TEST_PAR_SET_D2[i], obs[i]);
    }
}
