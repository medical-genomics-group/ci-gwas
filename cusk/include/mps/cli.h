#pragma once

#include <mps/marker_summary_stats.h>
#include <mps/marker_trait_summary_stats.h>
#include <mps/trait_summary_stats.h>
#include <mps/parent_set.h>

bool WRITE_FULL_CORRMATS = false;

struct CuskssArgs {
    bool merged;
    bool hetcor;
    bool trait_only;
    bool two_stage;
    bool time_indexed;
    float alpha;
    float pearson_sample_size;
    int max_level_one;
    int max_level_two;
    int depth;
    int block_ix;
    std::string block_path;
    std::string marker_ixs_path;
    std::string mxm_path;
    std::string mxp_path;
    std::string mxp_se_path;
    std::string pxp_path;
    std::string pxp_se_path;
    std::string time_index_path;
    std::string outdir;
};

struct CuskssSquareInputs {
    std::vector<float> correlations;
    std::vector<float> sample_sizes;
};

CuskssSquareInputs make_square_cuskss_inputs(
    const MarkerSummaryStats &mxm,
    const MarkerTraitSummaryStats &mxp,
    const TraitSummaryStats &pxp,
    const float pearson_sample_size,
    const bool heterogeneous_sample_sizes
);

ReducedGC run_cusk(ReducedGC gc, float threshold, int max_depth, int max_level, std::vector<int> &time_index_traits);

void cuskss(const CuskssArgs args);

void cusk(int argc, char *argv[]);

void prep_bed(int argc, char *argv[]);

void make_blocks(int argc, char *argv[]);