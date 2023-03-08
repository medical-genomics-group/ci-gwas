#include <gtest/gtest.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <test_data/bed_marker_test_set.h>

// Demonstrate some basic assertions.
TEST(Hello, BasicAssertions)
{
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

auto corr_matrix_size(size_t num_markers) -> size_t
{
    return num_markers * (num_markers - 1) / 2;
}

auto sparse_corr_matrix_size(
    size_t num_markers,
    size_t num_phen,
    size_t corr_width) -> size_t
{
    return (corr_width + num_phen) * (num_phen + num_markers);
}

TEST(AntiDiagSums, ExpectedReturnVals)
{
    float sums[11];
    memset(sums, 0.0, sizeof(sums));
    marker_corr_mat_antidiag_sums(BMT2_NUM_MARKERS, bmt2_marker_corrs, sums);
    for (size_t i = 0; i < 11; i++)
    {
        EXPECT_NEAR(sums[i], bmt2_marker_corr_antidiag_sums[i], 0.00001);
    }
}

TEST(SparseCorr, ExpectedReturnVals)
{
    const size_t corr_width = 3;
    const size_t batch_size = 5;
    const size_t cm_size = sparse_corr_matrix_size(BMT2_NUM_MARKERS, BMT2_NUM_PHEN, corr_width);

    float corrs[cm_size];
    memset(corrs, 0.0, sizeof(corrs));

    cu_corr_pearson_npn_batched_sparse(
        bmt2_marker_vals,
        bmt2_phen_vals,
        BMT2_NUM_MARKERS,
        BMT2_NUM_INDIVIDUALS,
        BMT2_NUM_PHEN,
        bmt2_marker_mean,
        bmt2_marker_std,
        corr_width,
        batch_size,
        corrs);

    for (size_t i = 0; i < cm_size; i++)
    {
        EXPECT_NEAR(corrs[i], bmt2_sparse_corrs[i], 0.00001);
    }
}

TEST(RowColIxFromLinear, ExpectedReturnVals3Markers)
{
    size_t row_ix, col_ix;
    size_t num_rows = 3;

    cu_ix_from_linear(0, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 0);
    EXPECT_EQ(col_ix, 1);

    cu_ix_from_linear(1, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 0);
    EXPECT_EQ(col_ix, 2);

    cu_ix_from_linear(2, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 1);
    EXPECT_EQ(col_ix, 2);
}

TEST(RowColIxFromLinear, ExpectedReturnVals7603Markers)
{
    size_t row_ix, col_ix;
    size_t num_rows = 7603;

    cu_ix_from_linear(28'899'002, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 7601);
    EXPECT_EQ(col_ix, 7602);

    cu_ix_from_linear(28'898'997, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 7599);
    EXPECT_EQ(col_ix, 7600);

    cu_ix_from_linear(28'898'996, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 7598);
    EXPECT_EQ(col_ix, 7602);

    cu_ix_from_linear(28'898'995, num_rows, &row_ix, &col_ix);
    EXPECT_EQ(row_ix, 7598);
    EXPECT_EQ(col_ix, 7601);
}

TEST(CuMarkerCorrPearsonBatched, ExpectedReturnVals)
{
    const size_t num_markers = BMT2_NUM_MARKERS;
    const size_t num_individuals = BMT2_NUM_INDIVIDUALS;
    const size_t marker_cm_size = corr_matrix_size(num_markers);
    const size_t stripe_width = 3;

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));

    cu_marker_corr_pearson_batched(
        bmt2_marker_vals,
        num_markers,
        num_individuals,
        bmt2_marker_mean,
        bmt2_marker_std,
        stripe_width,
        marker_corr);

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], bmt2_marker_corrs_pearson[i], 0.00001);
    }
}

TEST(CuMarkerCorrPearsonNpnBatched, ExpectedReturnVals)
{
    const size_t num_markers = BMT2_NUM_MARKERS;
    const size_t num_individuals = BMT2_NUM_INDIVIDUALS;
    const size_t marker_cm_size = corr_matrix_size(num_markers);
    const size_t stripe_width = 3;

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));

    cu_marker_corr_pearson_npn_batched(
        bmt2_marker_vals, num_markers, num_individuals, stripe_width, marker_corr);

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], bmt2_marker_corrs[i], 0.00001);
    }
}

TEST(CuCorrPearsonNpn, ExpectedReturnVals)
{
    const size_t num_markers = BMT_NUM_MARKERS;
    const size_t num_individuals = BMT_NUM_INDIVIDUALS;
    const size_t num_phen = BMT_NUM_PHEN;
    const size_t marker_cm_size = corr_matrix_size(num_markers);
    const size_t marker_phen_cm_size = num_markers * num_phen;
    const size_t phen_cm_size = corr_matrix_size(num_phen);

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));
    float marker_phen_corr[marker_phen_cm_size];
    memset(marker_phen_corr, 0.0, sizeof(marker_phen_corr));
    float phen_corr[phen_cm_size];
    memset(phen_corr, 0.0, sizeof(phen_corr));

    cu_corr_pearson_npn(
        bmt_marker_vals,
        bmt_phen_vals,
        num_markers,
        num_individuals,
        num_phen,
        bmt_marker_mean,
        bmt_marker_std,
        marker_corr,
        marker_phen_corr,
        phen_corr);

    float marker_corr_expected[marker_cm_size] = {0.299969, 0.924167, 0.0};
    float marker_phen_corr_expected[marker_phen_cm_size] = {
        -0.982847, -0.856312, -0.380668, -0.405905, -0.695899, -0.928034};
    float phen_corr_expected[phen_cm_size] = {0.808083};

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], marker_corr_expected[i], 0.00001);
    }

    for (size_t i = 0; i < marker_phen_cm_size; i++)
    {
        EXPECT_NEAR(marker_phen_corr[i], marker_phen_corr_expected[i], 0.00001);
    }

    for (size_t i = 0; i < phen_cm_size; i++)
    {
        EXPECT_NEAR(phen_corr[i], phen_corr_expected[i], 0.00001);
    }
}

TEST(CuCorrNpnBatched, ExpectedReturnVals)
{
    const size_t num_markers = BMT2_NUM_MARKERS;
    const size_t num_individuals = BMT2_NUM_INDIVIDUALS;
    const size_t num_phen = BMT2_NUM_PHEN;
    const size_t marker_cm_size = corr_matrix_size(num_markers);
    const size_t marker_phen_cm_size = num_markers * num_phen;
    const size_t phen_cm_size = corr_matrix_size(num_phen);
    const size_t stripe_width = 3;

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));
    float marker_phen_corr[marker_phen_cm_size];
    memset(marker_phen_corr, 0.0, sizeof(marker_phen_corr));
    float phen_corr[phen_cm_size];
    memset(phen_corr, 0.0, sizeof(phen_corr));

    cu_corr_pearson_npn_batched(
        bmt2_marker_vals,
        bmt2_phen_vals,
        num_markers,
        num_individuals,
        num_phen,
        bmt2_marker_mean,
        bmt2_marker_std,
        stripe_width,
        marker_corr,
        marker_phen_corr,
        phen_corr);

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], bmt2_marker_corrs[i], 0.00001);
    }

    for (size_t i = 0; i < marker_phen_cm_size; i++)
    {
        EXPECT_NEAR(marker_phen_corr[i], bmt2_marker_phen_corrs_pearson[i], 0.00001);
    }

    for (size_t i = 0; i < phen_cm_size; i++)
    {
        EXPECT_NEAR(phen_corr[i], bmt2_phen_corrs[i], 0.00001);
    }
}

TEST(CuMarkerPearson, ExpectedReturnVals3Markers)
{
    const size_t num_markers = BMT_NUM_MARKERS;
    const size_t num_individuals = BMT_NUM_INDIVIDUALS;
    const size_t marker_cm_size = corr_matrix_size(num_markers);

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));

    cu_marker_corr_pearson(
        bmt_marker_vals, num_markers, num_individuals, bmt_marker_mean, bmt_marker_std,
        marker_corr);

    float marker_corr_expected[marker_cm_size] = {0.2540839, 0.80403025, 0.04012862};

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], marker_corr_expected[i], 0.00001);
    }
}

TEST(CuMarkerPearson, ExpectedReturnVals7Markers)
{
    const size_t num_markers = BMT2_NUM_MARKERS;
    const size_t num_individuals = BMT2_NUM_INDIVIDUALS;
    const size_t marker_cm_size = corr_matrix_size(num_markers);

    float marker_corr[marker_cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));

    cu_marker_corr_pearson(
        bmt2_marker_vals,
        num_markers,
        num_individuals,
        bmt2_marker_mean,
        bmt2_marker_std,
        marker_corr);

    for (size_t i = 0; i < marker_cm_size; i++)
    {
        EXPECT_NEAR(marker_corr[i], bmt2_marker_corrs_pearson[i], 0.00001);
    }
}
