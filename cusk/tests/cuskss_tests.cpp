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

void expect_eq_vec(std::vector<std::string> exp, std::vector<std::string> obs)
{
    EXPECT_EQ(exp.size(), obs.size());
    for (size_t i = 0; i < exp.size(); i++)
    {
        EXPECT_EQ(exp[i], obs[i]);
    }
}

TEST(cuskss_pearson, expected_results)
{
    std::cout << "TESTPRINT" << std::endl;
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
    
    std::vector<int> exp_adj = {
        1, 1, 1,
        1, 1, 1,
        1, 1, 1
    };
    std::vector<int> obs_adj = read_ints_from_binary(temp_dir / "trait_only.adj");

    expect_eq_vec(exp_adj, obs_adj);
    // "/tmp/cuskss_test_files/trait_only.corr"
    // "/tmp/cuskss_test_files/trait_only.ixs"
    // "/tmp/cuskss_test_files/trait_only.adj"
    // "/tmp/cuskss_test_files/trait_only.mdim"
    // for (const auto & entry : std::filesystem::directory_iterator(temp_dir))
    //     std::cout << entry.path() << std::endl;
    // which files do we expect?
    // what file contents do we expect?
    std::filesystem::remove_all(temp_dir);
}