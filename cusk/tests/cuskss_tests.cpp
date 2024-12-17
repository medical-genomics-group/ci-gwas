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

TEST(cuskss_pearson, expected_results)
{
    std::cout << "TESTPRINT" << std::endl;
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "cuskss_test_files";
    std::filesystem::create_directories(temp_dir);
    std::string marker_ixs_path = "../../tests/test_files/marker_indices.bin";
    std::string mxm_path = "../../tests/test_files/small_mxm.bin";
    std::string mxp_path = "../../tests/test_files/marker_trait_summary_stats.txt";
    std::string pxp_path = "../../tests/test_files/trait_summary_stats.txt";
    CuskssArgs args = {
        false,              // merged
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
    for (const auto & entry : std::filesystem::directory_iterator(temp_dir))
        std::cout << entry.path() << std::endl;
    // which files do we expect?
    // what file contents do we expect?
    std::filesystem::remove(temp_dir);
}