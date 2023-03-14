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

TEST(PrepBim, ExpectedReturnVals)
{
    bimInfo exp;
    exp.number_of_lines = 5;
    exp.chr_ids = {"1", "19"};
    exp.num_markers_on_chr = {3, 2};

    bimInfo obs = prep_bim("../../tests/test_files/small.bim", "../../tests/test_files");

    EXPECT_EQ(obs.number_of_lines, exp.number_of_lines);

    for (size_t i = 0; i < 2; ++i)
    {
        EXPECT_EQ(obs.chr_ids[i], exp.chr_ids[i]);
        EXPECT_EQ(obs.num_markers_on_chr[i], exp.num_markers_on_chr[i]);
    }

    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.bim"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.bim"));
}

TEST(ParseBedDeath, WrongMagicNumberOne)
{
    ASSERT_DEATH(
        {
            prep_bed(
                "../../tests/test_files/wrong_magic_num_one.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        },
        "unexpected magic number in bed file.");
}

TEST(ParseBedDeath, WrongMagicNumberTwo)
{
    ASSERT_DEATH(
        {
            prep_bed(
                "../../tests/test_files/wrong_magic_num_two.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        },
        "unexpected magic number in bed file.");
}

TEST(ParseBedDeath, WrongMagicNumberThree)
{
    ASSERT_DEATH(
        {
            prep_bed(
                "../../tests/test_files/wrong_magic_num_three.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        },
        "unexpected magic number in bed file.");
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
    prep_bed(
        "../../tests/test_files/small.bed",
        "../../tests/test_files/small.bim",
        "../../tests/test_files/small.fam",
        "../../tests/test_files",
        1);

    // .bed
    std::vector<unsigned char> one_bed_content = read_file("../../tests/test_files/1.bed");
    std::vector<unsigned char> one_bed_exp = {0xea, 0x8f, 0x0f, 0x38, 0x8e, 0x03, 0xea, 0x8a, 0x0f};
    EXPECT_EQ(one_bed_content.size(), one_bed_exp.size());
    for (size_t i = 0; i < one_bed_exp.size(); ++i)
    {
        EXPECT_EQ(one_bed_exp[i], one_bed_content[i]);
    }

    std::vector<unsigned char> nt_bed_content = read_file("../../tests/test_files/19.bed");
    std::vector<unsigned char> nt_bed_exp = {0xb2, 0x2c, 0x0b, 0xcb, 0xb2, 0x0c};
    EXPECT_EQ(nt_bed_content.size(), nt_bed_exp.size());
    for (size_t i = 0; i < nt_bed_exp.size(); ++i)
    {
        EXPECT_EQ(nt_bed_exp[i], nt_bed_content[i]);
    }

    // .stds
    std::vector<float> one_stds_exp = {0.66332496, 0.83066239, 0.6};
    std::vector<float> one_stds_obs = read_floats_from_lines("../../tests/test_files/1.stds");
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(one_stds_exp[i], one_stds_obs[i], 0.000001);
    }

    std::vector<float> nt_stds_exp = {0.77459667, 0.83066239};
    std::vector<float> nt_stds_obs = read_floats_from_lines("../../tests/test_files/19.stds");
    for (size_t i = 0; i < 2; ++i)
    {
        EXPECT_NEAR(nt_stds_exp[i], nt_stds_obs[i], 0.000001);
    }

    // .means
    std::vector<float> one_means_exp = {0.6, 1.1, 0.8};
    std::vector<float> one_means_obs = read_floats_from_lines("../../tests/test_files/1.means");
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(one_means_exp[i], one_means_obs[i], 0.000001);
    }

    std::vector<float> nt_means_exp = {1., 0.9};
    std::vector<float> nt_means_obs = read_floats_from_lines("../../tests/test_files/19.means");
    for (size_t i = 0; i < 2; ++i)
    {
        EXPECT_NEAR(nt_means_exp[i], nt_means_obs[i], 0.000001);
    }

    // .dim
    std::vector<int> one_dims_exp = {10, 3};
    std::vector<int> one_dims_obs = read_ints_from_lines("../../tests/test_files/1.dim");
    for (size_t i = 0; i < 2; ++i)
    {
        EXPECT_EQ(one_dims_exp[i], one_dims_obs[i]);
    }

    std::vector<int> nt_dims_exp = {10, 2};
    std::vector<int> nt_dims_obs = read_ints_from_lines("../../tests/test_files/19.dim");
    for (size_t i = 0; i < 2; ++i)
    {
        EXPECT_EQ(nt_dims_exp[i], nt_dims_obs[i]);
    }

    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.dim"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.bed"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.means"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.stds"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.bed"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.dim"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.means"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.stds"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.bim"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.bim"));
}
