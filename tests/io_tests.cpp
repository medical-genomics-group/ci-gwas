#include <gtest/gtest.h>
#include <mps/io.h>
#include <vector>
#include <string>

TEST(LoadBinaryTest, ExpectedResults) {
    size_t nbytes = 5;
    std::vector<unsigned char> exp = {0x6c, 0x1b, 0x00, 0xea, 0x8f};
    std::vector<unsigned char> res(nbytes);
    read_n_bytes_from_binary("../../tests/test_files/small.bed", nbytes, res);
    EXPECT_EQ(res.size(), nbytes);
    for (size_t i = 0; i < nbytes; ++i) {
        EXPECT_EQ(exp[i], res[i]);
    }
}
