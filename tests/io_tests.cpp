#include <gtest/gtest.h>
#include <mps/io.h>
#include <vector>
#include <string>
#include <filesystem>

TEST(LoadBinaryTest, ExpectedResults) {
    size_t nbytes = 5;
    std::vector<unsigned char> exp = {0x6c, 0x1b, 0x01, 0xea, 0x8f};
    std::vector<unsigned char> res(nbytes);
    read_n_bytes_from_binary("../../tests/test_files/small.bed", nbytes, res);
    
    EXPECT_EQ(res.size(), nbytes);
    for (size_t i = 0; i < nbytes; ++i) {
        EXPECT_EQ(exp[i], res[i]);
    }
}

TEST(ReadFloatsFromLinesTest, ExpectedResults) {
    std::vector<float> exp_vals = {0.72031609, 0.59822862, 0.47614114, -0.62264611,
                                   -1.110996, -1.23308348, 2.06327829, 0.59822862,
                                   -0.86682106, -0.62264611};
    std::vector<float> read_vals = {};
    std::string loc = "../../tests/test_files/small.phen";
    read_floats_from_lines(loc, read_vals);
    
    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(exp_vals[i], read_vals[i]);
    }
}

TEST(FloatToAndFromBinaryTest, ExpectedResults) {
    std::vector<float> data = {0.0, -1.1, 2.2, -3.3, 4.4};
    std::string loc = "../../tests/test_files/floats.bin";
    write_floats_to_binary(data.data(), 4, loc);
    std::vector<float> read = read_floats_from_binary(loc);
    
    EXPECT_EQ(read.size(), 4);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(read[i], data[i]);
    }

    EXPECT_TRUE(std::filesystem::remove(loc));
}
