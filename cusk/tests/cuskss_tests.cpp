#include <gtest/gtest.h>
#include <mps/bim.h>
#include <mps/io.h>
#include <mps/marker_summary_stats.h>
#include <mps/marker_trait_summary_stats.h>
#include <mps/trait_summary_stats.h>
#include <mps/cli.h>

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

void expect_eq_vec(std::vector<int> exp, std::vector<int> obs)
{
    EXPECT_EQ(exp.size(), obs.size());
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(cuskss, trait_only_merged_expected_results)
{
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "cuskss_test_files";
    std::filesystem::create_directories(temp_dir);
    std::string marker_ixs_path = "../../tests/test_files/marker_indices.bin";
    std::string mxm_path = "../../tests/test_files/small_mxm.bin";
    std::string mxp_path = "../../tests/test_files/marker_trait_summary_stats.txt";
    std::string pxp_path = "../../tests/test_files/trait_summary_stats.txt";
    check_path(pxp_path);
    CuskssArgs args = {
        true,               // merged
        false,              // hetcor
        true,               // trait_only
        false,              // two_stage
        false,              // time_indexed
        0.0001,             // alpha
        500000.f,           // pearson_sample_size
        3,                  // max_level
        0,                  // max_level_two
        1,                  // depth
        0,                  // block_ix
        "NULL",             // block_path
        marker_ixs_path,    // marker_ixs_path
        mxm_path,           // mxm_path
        mxp_path,           // mxp_path
        "NULL",             // mxp_se_path
        pxp_path,           // pxp_path
        "NULL",             // pxp_se_path
        "NULL",             // time_index_path
        temp_dir            // outdir
    };
    cuskss(args);
    
    // output adjacency matrix as expected
    std::vector<int> exp_adj = {
        0, 1, 1,
        1, 0, 1,
        1, 1, 0
    };
    std::vector<int> obs_adj = read_ints_from_binary(temp_dir / "trait_only.adj");
    expect_eq_vec(exp_adj, obs_adj);

    // output adjacency matrix as expected
    std::vector<float> exp_corr = {
        1.0, 0.0608594558771734, 0.074239793758568,
        0.0608594558771734, 1.0, 0.0675875270156859,
        0.074239793758568, 0.0675875270156859, 1.0
    };
    std::vector<float> obs_corr = read_floats_from_binary(temp_dir / "trait_only.corr");
    expect_near_vec(exp_corr, obs_corr);
    
    std::filesystem::remove_all(temp_dir);
}

TEST(cuskss, pearson_two_stage_merged_expected_results)
{
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "cuskss_test_files";
    std::filesystem::create_directories(temp_dir);
    std::string marker_ixs_path = "../../tests/test_files/marker_indices.bin";
    std::string mxm_path = "../../tests/test_files/small_mxm.bin";
    std::string mxp_path = "../../tests/test_files/marker_trait_summary_stats.txt";
    std::string pxp_path = "../../tests/test_files/trait_summary_stats.txt";
    check_path(pxp_path);
    CuskssArgs args = {
        true,               // merged
        false,              // hetcor
        false,              // trait_only
        true,               // two_stage
        false,              // time_indexed
        0.0001,             // alpha
        500000.f,           // pearson_sample_size
        3,                  // max_level
        1,                  // max_level_two
        1,                  // depth
        0,                  // block_ix
        "NULL",             // block_path
        marker_ixs_path,    // marker_ixs_path
        mxm_path,           // mxm_path
        mxp_path,           // mxp_path
        "NULL",             // mxp_se_path
        pxp_path,           // pxp_path
        "NULL",             // pxp_se_path
        "NULL",             // time_index_path
        temp_dir            // outdir
    };
    cuskss(args);
        
   // selected markers as expected
    std::vector<int> exp_ixs = {2, 3, 4, 5};
    std::vector<int> obs_ixs = read_ints_from_binary(temp_dir / "cuskss_merged.ixs");
    expect_eq_vec(exp_ixs, obs_ixs);
    
    // output adjacency matrix as expected
    std::vector<int> exp_adj = {
        0, 0, 0, 1,
        0, 0, 1, 1,
        0, 1, 0, 1,
        1, 1, 1, 0,
    };
    std::vector<int> obs_adj = read_ints_from_binary(temp_dir / "cuskss_merged.adj");
    expect_eq_vec(exp_adj, obs_adj);

    // output adjacency matrix as expected
    std::vector<float> exp_corr = {
        1.0, 0.0005, 0.0001, -0.01,
        0.0005, 1.0, 0.0608594558771734, 0.074239793758568,
        0.0001, 0.0608594558771734, 1.0, 0.0675875270156859,
        -0.01, 0.074239793758568, 0.0675875270156859, 1.0,
        
    };
    std::vector<float> obs_corr = read_floats_from_binary(temp_dir / "cuskss_merged.corr");
    expect_near_vec(exp_corr, obs_corr);
    std::filesystem::remove_all(temp_dir);
}

TEST(cuskss, pearson_two_stage_block_expected_results)
{
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "cuskss_test_files";
    std::filesystem::create_directories(temp_dir);
    std::string marker_ixs_path = "../../tests/test_files/marker_indices.bin";
    std::string mxm_path = "../../tests/test_files/small_mxm.bin";
    std::string mxp_path = "../../tests/test_files/marker_trait_summary_stats.txt";
    std::string pxp_path = "../../tests/test_files/trait_summary_stats.txt";
    std::string block_path = "../../tests/test_files/blocks.txt";
    std::string time_index_path = "../../tests/test_files/time_index.txt";
    check_path(pxp_path);
    CuskssArgs args = {
        false,              // merged
        false,              // hetcor
        false,              // trait_only
        true,               // two_stage
        false,              // time_indexed
        0.0001,             // alpha
        500000.f,           // pearson_sample_size
        3,                  // max_level
        1,                  // max_level_two
        1,                  // depth
        1,                  // block_ix
        block_path,         // block_path
        marker_ixs_path,    // marker_ixs_path
        mxm_path,           // mxm_path
        mxp_path,           // mxp_path
        "NULL",             // mxp_se_path
        pxp_path,           // pxp_path
        "NULL",             // pxp_se_path
        time_index_path,    // time_index_path
        temp_dir            // outdir
    };
    cuskss(args);
        
   // selected markers as expected
    std::vector<int> exp_ixs = {2, 3, 4, 5};
    std::vector<int> obs_ixs = read_ints_from_binary(temp_dir / "1_2_2.ixs");
    expect_eq_vec(exp_ixs, obs_ixs);
    
    // output adjacency matrix as expected
    std::vector<int> exp_adj = {
        0, 0, 0, 1,
        0, 0, 1, 1,
        0, 1, 0, 1,
        1, 1, 1, 0,
    };
    std::vector<int> obs_adj = read_ints_from_binary(temp_dir / "1_2_2.adj");
    expect_eq_vec(exp_adj, obs_adj);

    // output adjacency matrix as expected
    std::vector<float> exp_corr = {
        1.0, 0.0005, 0.0001, -0.01,
        0.0005, 1.0, 0.0608594558771734, 0.074239793758568,
        0.0001, 0.0608594558771734, 1.0, 0.0675875270156859,
        -0.01, 0.074239793758568, 0.0675875270156859, 1.0,
        
    };
    std::vector<float> obs_corr = read_floats_from_binary(temp_dir / "1_2_2.corr");
    expect_near_vec(exp_corr, obs_corr);
    std::filesystem::remove_all(temp_dir);
}