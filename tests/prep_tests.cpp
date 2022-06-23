#include <gtest/gtest.h>
#include <mps/prep_markers.h>


TEST(CountLinesTest, ExpectedReturnVals) {
    EXPECT_EQ(count_lines("../../tests/test_files/small.fam"), 10);
}

