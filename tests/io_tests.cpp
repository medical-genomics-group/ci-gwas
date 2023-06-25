#include <gtest/gtest.h>
#include <mps/bim.h>
#include <mps/io.h>
#include <mps/marker_summary_stats.h>
#include <mps/marker_trait_summary_stats.h>
#include <mps/trait_summary_stats.h>

#include <filesystem>
#include <string>
#include <vector>

void expect_near_vec(std::vector<float> exp, std::vector<float> obs)
{
    EXPECT_EQ(exp.size(), obs.size());
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_NEAR(exp[i], obs[i], 0.001);
    }
}

void expect_eq_vec(std::vector<std::string> exp, std::vector<std::string> obs)
{
    EXPECT_EQ(exp.size(), obs.size());
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(read_trait_summary_stats, expected_results)
{
    std::vector<std::string> exp_header = {"AT", "BMI", "CAD"};

    std::vector<float> exp_corrs = {
        0.0,
        0.0608594558771734,
        0.074239793758568,
        0.0608594558771734,
        0.0,
        0.0675875270156859,
        0.074239793758568,
        0.0675875270156859,
        0.0};

    int exp_num_phen = 3;

    TraitSummaryStats obs = TraitSummaryStats("../../tests/test_files/trait_summary_stats.txt");

    EXPECT_EQ(exp_num_phen, obs.num_phen);
    expect_eq_vec(exp_num_phen, obs.header);
    expect_near_vec(exp_header, obs.corrs);
}

TEST(read_correlations_from_mtx, expected_results)
{
    std::vector<float> exp_vals = {
        1.0,
        .01676169840066658,
        -.026459418125921043,
        .01676169840066658,
        1.0,
        .013313834590098538,
        -.026459418125921043,
        .013313834590098538,
        1.0};

    std::vector<float> obs_vals = read_correlations_from_mtx("../../tests/test_files/test.mtx");

    EXPECT_EQ(obs_vals.size(), exp_vals.size());
    for (size_t i = 0; i < exp_vals.size(); ++i)
    {
        EXPECT_EQ(exp_vals[i], obs_vals[i]);
    }
}

TEST(LoadBinary, ExpectedResults)
{
    size_t nbytes = 5;
    std::vector<unsigned char> exp = {0x6c, 0x1b, 0x01, 0xea, 0x8f};
    std::vector<unsigned char> res(nbytes);
    read_n_bytes_from_binary("../../tests/test_files/small.bed", nbytes, res);

    EXPECT_EQ(res.size(), nbytes);
    for (size_t i = 0; i < nbytes; ++i)
    {
        EXPECT_EQ(exp[i], res[i]);
    }
}

TEST(ReadFloatsFromLines, ExpectedResults)
{
    std::vector<float> exp_vals = {
        0.72031609,
        0.59822862,
        0.47614114,
        -0.62264611,
        -1.110996,
        -1.23308348,
        2.06327829,
        0.59822862,
        -0.86682106,
        -0.62264611};
    std::vector<float> read_vals = {};
    std::string loc = "../../tests/test_files/small.phen";
    read_floats_from_lines(loc, read_vals);

    for (size_t i = 0; i < 10; ++i)
    {
        EXPECT_EQ(exp_vals[i], read_vals[i]);
    }
}

TEST(FloatToAndFromBinary, ExpectedResults)
{
    std::vector<float> data = {0.0, -1.1, 2.2, -3.3, 4.4};
    std::string loc = "../../tests/test_files/floats.bin";
    write_floats_to_binary(data.data(), 4, loc);
    std::vector<float> read = read_floats_from_binary(loc);

    EXPECT_EQ(read.size(), 4);
    for (size_t i = 0; i < 4; ++i)
    {
        EXPECT_EQ(read[i], data[i]);
    }

    EXPECT_TRUE(std::filesystem::remove(loc));
}

TEST(LoadBlocks, ExpectedResults)
{
    std::vector<MarkerBlock> exp{
        MarkerBlock("1", 101, 202), MarkerBlock("1", 303, 404), MarkerBlock("1", 505, 606)};
    std::vector<MarkerBlock> obs = read_blocks_from_file("../../tests/test_files/test.blocks");
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(LoadBedBlock, ExpectedResults)
{
    MarkerBlock block("1", 1, 2);
    BedDims dims(10, 5);
    BimInfo bim("../../tests/test_files/small.bim");
    std::vector<unsigned char> obs =
        read_block_from_bed("../../tests/test_files/small.bed", block, dims, bim);
    std::vector<unsigned char> exp{0x38, 0x8e, 0xf3, 0xea, 0x8a, 0xff};
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(LoadBedBlockSecondChr, ExpectedResults)
{
    MarkerBlock block("19", 0, 1);
    BedDims dims(10, 5);
    BimInfo bim("../../tests/test_files/small.bim");
    std::vector<unsigned char> obs =
        read_block_from_bed("../../tests/test_files/small.bed", block, dims, bim);
    std::vector<unsigned char> exp{0xb2, 0x2c, 0xfb, 0xcb, 0xb2, 0xfc};
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}
