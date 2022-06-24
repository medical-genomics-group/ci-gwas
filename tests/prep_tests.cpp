#include <gtest/gtest.h>
#include <mps/prep_markers.h>
#include <string>

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
