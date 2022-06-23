#include <gtest/gtest.h>
#include <mps/prep_markers.h>
#include <string>

TEST(CountLinesTest, ExpectedReturnVals) {
    EXPECT_EQ(count_lines("../../tests/test_files/small.fam"), 10);
}

TEST(SplitLineTest, ExpectedReturnVals) {
    std::string bim_line[6];
    std::string str = "19      rs150203985     0       240963  A       G";
    std::string exp_fields[6] = {"19","rs150203985", "0", "240963", "A", "G"};
    split_line(str, bim_line);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(bim_line[i], exp_fields[i]);
    }
}

//TEST(ParseBimTest, ExpectedReturnVals) {
//}
