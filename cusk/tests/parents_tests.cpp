#include <gtest/gtest.h>
#include <mps/parent_set.h>
#include <test_data/parent_set_test_set.h>

TEST(subset_variables, ExpectedResultsDepth0)
{
    std::unordered_set<int> obs_set =
        subset_variables(TEST_ADJ_MAT, TEST_NUM_MARKERS + TEST_NUM_PHEN, TEST_NUM_MARKERS, 0);
    std::vector<int> obs = set_to_vec(obs_set);
    ASSERT_EQ(obs.size(), TEST_PAR_SET_D0.size());
    for (size_t i = 0; i < obs.size(); ++i)
    {
        EXPECT_EQ(TEST_PAR_SET_D0[i], obs[i]);
    }
}

TEST(subset_variables, ExpectedResultsDepth1)
{
    std::unordered_set<int> obs_set =
        subset_variables(TEST_ADJ_MAT, TEST_NUM_MARKERS + TEST_NUM_PHEN, TEST_NUM_MARKERS, 1);
    std::vector<int> obs = set_to_vec(obs_set);
    ASSERT_EQ(obs.size(), TEST_PAR_SET_D1.size());
    for (size_t i = 0; i < obs.size(); ++i)
    {
        EXPECT_EQ(TEST_PAR_SET_D1[i], obs[i]);
    }
}

TEST(subset_variables, ExpectedResultsDepth2)
{
    std::unordered_set<int> obs_set =
        subset_variables(TEST_ADJ_MAT, TEST_NUM_MARKERS + TEST_NUM_PHEN, TEST_NUM_MARKERS, 2);
    std::vector<int> obs = set_to_vec(obs_set);
    ASSERT_EQ(obs.size(), TEST_PAR_SET_D2.size());
    for (size_t i = 0; i < obs.size(); ++i)
    {
        EXPECT_EQ(TEST_PAR_SET_D2[i], obs[i]);
    }
}
