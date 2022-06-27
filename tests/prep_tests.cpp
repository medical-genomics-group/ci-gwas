#include <gtest/gtest.h>
#include <mps/prep_markers.h>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>

TEST(CountLinesTest, ExpectedReturnVals) {
    EXPECT_EQ(count_lines("../../tests/test_files/small.fam"), 10);
}

TEST(SplitBimLineTest, ExpectedReturnVals) {
    std::string bim_line[6];
    std::string str = "19      rs150203985     0       240963  A       G";
    std::string exp_fields[6] = {"19","rs150203985", "0", "240963", "A", "G"};
    split_bim_line(str, bim_line);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(bim_line[i], exp_fields[i]);
    }
}

TEST(ParseBimTest, ExpectedReturnVals) {
    bimInfo exp;
    exp.number_of_lines = 5;
    exp.chr_ids = {"1", "19"};
    exp.num_markers_on_chr = {3, 2};

    bimInfo obs = parse_bim("../../tests/test_files/small.bim");
   
    EXPECT_EQ(obs.number_of_lines, exp.number_of_lines);

    for (size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(obs.chr_ids[i], exp.chr_ids[i]);
        EXPECT_EQ(obs.num_markers_on_chr[i], exp.num_markers_on_chr[i]);
    }
}

TEST(ParseBedDeathTest, WrongMagicNumberOne) {
    ASSERT_DEATH({
        prep_bed(
                "../../tests/test_files/wrong_magic_num_one.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        }, "unexpected magic number in bed file.");
}

TEST(ParseBedDeathTest, WrongMagicNumberTwo) {
    ASSERT_DEATH({
        prep_bed(
                "../../tests/test_files/wrong_magic_num_two.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        }, "unexpected magic number in bed file.");
}

TEST(ParseBedDeathTest, WrongMagicNumberThree) {
    ASSERT_DEATH({
        prep_bed(
                "../../tests/test_files/wrong_magic_num_three.bed",
                "../../tests/test_files/small.bim",
                "../../tests/test_files/small.fam",
                "../../tests/test_files",
                1);
        }, "unexpected magic number in bed file.");
}

auto read_file(const std::string filename) -> std::vector<unsigned char>
{
    // open the file:
    std::ifstream file(filename, std::ios::binary);

    // read the data:
    return std::vector<unsigned char>(
            (std::istreambuf_iterator<char>(file)),
            std::istreambuf_iterator<char>());
}

TEST(ParseBed, CorrectOutFilesGenerated) {
    prep_bed(
        "../../tests/test_files/small.bed",
        "../../tests/test_files/small.bim",
        "../../tests/test_files/small.fam",
        "../../tests/test_files",
        1);

    std::vector<unsigned char> one_bed_content = read_file("../../tests/test_files/1.bed");
    std::vector<unsigned char> one_bed_exp = {0xea, 0x8f, 0x0f, 0x38, 0x8e, 0x03, 0xea, 0x8a, 0x0f};
    EXPECT_EQ(one_bed_content.size(), one_bed_exp.size());
    for (size_t i = 0; i < one_bed_exp.size(); ++i) {
        EXPECT_EQ(one_bed_exp[i], one_bed_content[i]);
    }

    std::vector<unsigned char> nt_bed_content = read_file("../../tests/test_files/19.bed");
    std::vector<unsigned char> nt_bed_exp = {0xb2, 0x2c, 0x0b, 0xcb, 0xb2, 0x0c};
    EXPECT_EQ(nt_bed_content.size(), nt_bed_exp.size());
    for (size_t i = 0; i < nt_bed_exp.size(); ++i) {
        EXPECT_EQ(nt_bed_exp[i], nt_bed_content[i]);
    }

    

    // bed 1 should contain \xea, \x8f, \xff, \x38, \x8e, \xf3, \xea, \x8a, \xff
    // bed 19 should contain \xb2, \x2c, \xfb, \xcb, \xb2, \xfc
    // std should be [0.66332496, 0.83066239, 0.6       , 0.77459667, 0.83066239]
    // mean shuold be [0.6, 1.1, 0.8, 1. , 0.9]

    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.bed"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.means"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/19.stds"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.bed"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.means"));
    EXPECT_TRUE(std::filesystem::remove("../../tests/test_files/1.stds"));
}
