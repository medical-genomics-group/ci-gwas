#include <gtest/gtest.h>
#include <mps/bim.h>
#include <mps/io.h>
#include <mps/prep.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

TEST(BimStat, MarkersInDistance)
{
    EXPECT_EQ(num_markers_within_distance("../../tests/test_files/distance.bim", 3000), 2);
}

TEST(CountLines, ExpectedReturnVals)
{
    EXPECT_EQ(count_lines("../../tests/test_files/small.fam"), 10);
}

TEST(SplitBimLine, ExpectedReturnVals)
{
    std::string bim_line[6];
    std::string str = "19      rs150203985     0       240963  A       G";
    std::string exp_fields[6] = {"19", "rs150203985", "0", "240963", "A", "G"};
    split_bim_line(str, bim_line);
    for (size_t i = 0; i < 6; ++i)
    {
        EXPECT_EQ(bim_line[i], exp_fields[i]);
    }
}

auto read_file(const std::string filename) -> std::vector<unsigned char>
{
    // open the file:
    std::ifstream file(filename, std::ios::binary);

    // read the data:
    return std::vector<unsigned char>(
        (std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
}

TEST(ParseBed, CorrectOutFilesGenerated)
{
    prep_bed_no_impute("../../tests/test_files/small");

    // .stds
    std::vector<float> one_stds_exp = {0.66332496, 0.83066239, 0.6, 0.77459667, 0.83066239};
    std::vector<float> one_stds_obs = read_floats_from_lines("../../tests/test_files/small.stds");
    for (size_t i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(one_stds_exp[i], one_stds_obs[i], 0.000001);
    }

    // .means
    std::vector<float> one_means_exp = {0.6, 1.1, 0.8, 1., 0.9};
    std::vector<float> one_means_obs = read_floats_from_lines("../../tests/test_files/small.means");
    for (size_t i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(one_means_exp[i], one_means_obs[i], 0.000001);
    }

    // .dim
    BedDims exp(10, 5);
    BedDims obs("../../tests/test_files/small.dim");
    EXPECT_EQ(exp, obs);

    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/small.dim"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/small.means"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/small.stds"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/small.modes"));
}
