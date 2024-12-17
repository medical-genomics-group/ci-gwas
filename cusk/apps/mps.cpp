#include <math.h>
#include <mps/bim.h>
#include <mps/blocking.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/cuPC-S.h>
#include <mps/hetcor-cuPC-S.h>
#include <mps/cuPC_call_prep.h>
#include <mps/io.h>
#include <mps/marker_summary_stats.h>
#include <mps/marker_trait_summary_stats.h>
#include <mps/parent_set.h>
#include <mps/phen.h>
#include <mps/prep.h>
#include <mps/trait_summary_stats.h>
#include <sys/stat.h>

#include <numeric>
#include <array>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

bool WRITE_FULL_CORRMATS = false;

ReducedGC run_cusk(
    ReducedGC gc,
    float threshold,
    int max_depth,
    int max_level, 
    std::vector<int> &time_index_traits
)
{
    std::vector<int> time_index_gc(gc.num_var, 0);
    int read_ix = 0;
    for (int i = gc.num_markers(); i < gc.num_var; i++) {
        time_index_gc[i] = time_index_traits[read_ix];
        read_ix++;
    }
    int num_var = gc.num_var;
    int start_level = 0;
    hetcor_skeleton(
        gc.C.data(),
        &num_var,
        gc.G.data(),
        gc.S.data(),
        &threshold,
        &start_level,
        &max_level,
        time_index_gc.data()
    );
    std::unordered_set<int> variable_subset =
        subset_variables(gc.G, gc.num_var, gc.num_markers(), max_depth);
    return reduce_gc(
        gc.G, gc.C, gc.S, variable_subset, gc.num_var, gc.num_phen, ML, gc.new_to_old_indices
    );
}

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
){
    // make n2 matrix to please cuPC
    size_t num_phen = pxp.get_num_phen();
    size_t num_markers = mxm.get_num_markers();
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs(num_var * num_var, 1.0);
    std::vector<float> sq_ess(num_var * num_var, pearson_sample_size);
    std::vector<float> marker_corr = mxm.get_corrs();

    size_t sq_row_ix = 0;
    size_t sq_col_ix = 0;
    for (size_t i = 0; i < marker_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)marker_corr[i];
        if (sq_col_ix == num_markers - 1)
        {
            ++sq_row_ix;
            // marker_corr is square
            sq_col_ix = 0;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    std::vector<float> marker_phen_corr = mxp.get_corrs();
    std::vector<float> marker_phen_sample_sizes;
    if (heterogeneous_sample_sizes) {
        std::vector<float> marker_phen_sample_sizes = mxp.get_sample_sizes();
    }
    sq_row_ix = 0;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < marker_phen_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)marker_phen_corr[i];
        sq_corrs[num_var * sq_col_ix + sq_row_ix] = (float)marker_phen_corr[i];
        if (heterogeneous_sample_sizes) {
            sq_ess[num_var * sq_row_ix + sq_col_ix] = (float)marker_phen_sample_sizes[i];
            sq_ess[num_var * sq_col_ix + sq_row_ix] = (float)marker_phen_sample_sizes[i];
        }
        if (sq_col_ix == (num_var - 1))
        {
            sq_col_ix = num_markers;
            ++sq_row_ix;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    std::vector<float> phen_sample_sizes;
    if (heterogeneous_sample_sizes) {
        std::vector<float> phen_sample_sizes = pxp.get_sample_sizes();
    }
    std::vector<float> phen_corr = pxp.get_corrs();
    sq_row_ix = num_markers;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < phen_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)phen_corr[i];
        if (heterogeneous_sample_sizes) {
            sq_ess[num_var * sq_row_ix + sq_col_ix] = (float)phen_sample_sizes[i];
        }
        if (sq_col_ix == (num_var - 1))
        {
            ++sq_row_ix;
            sq_col_ix = num_markers;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    CuskssSquareInputs result = {sq_corrs, sq_ess};
    return result;
}

/**
 * @brief Compute number of variables from upper triangular matrix (diag excluded)
 *
 * @return int
 */
auto num_variables_from_matrix_size(const size_t num_matrix_entries) -> size_t
{
    return (1 + (size_t)sqrt(1.0 + 8.0 * (long double)num_matrix_entries)) / 2;
}

void check_nargs(const int argc, const int nargs, const std::string usage)
{
    if (argc < nargs)
    {
        std::cout << usage << std::endl;
        exit(1);
    }
}

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

void cuskss(const CuskssArgs args)
{
    // load input data
    std::cout << "Loading input files" << std::endl;
    
    if (args.merged) {
        std::cout << "Loading marker indices" << std::endl;
        std::vector<int> marker_ixs = read_ints_from_binary(args.marker_ixs_path);
    } else {
        std::cout << "Loading block file" << std::endl;
        std::vector<MarkerBlock> blocks = read_blocks_from_file(args.block_path);
        MarkerBlock block = blocks[args.block_ix];
    }

    TraitSummaryStats pxp;
    std::cout << "Loading pxp" << std::endl;
    if (args.hetcor) {
        TraitSummaryStats pxp = TraitSummaryStats(args.pxp_path, args.pxp_se_path);
    } else {
        TraitSummaryStats pxp = TraitSummaryStats(args.pxp_path, args.pearson_sample_size);
    }
    size_t num_phen = pxp.get_num_phen();

    if (args.time_indexed) {
        std::cout << "Loading time_indices" << std::endl;
        std::vector<int> time_index_traits = read_ints_from_lines(args.time_index_path);
    } else {
        // there is no time index file to load, just put all traits at time 1
        std::vector<int> time_index_traits(num_phen, 1);
    }
    
    if (args.trait_only) {
        size_t num_markers = 0;
        size_t num_var = num_markers + num_phen;
        std::vector<float> sq_corrs = pxp.get_corrs();
        std::vector<float> sq_ess = pxp.get_sample_sizes();

        // init new-to-old indices for initial gc
        std::vector<int> nto_ixs(num_var);
        std::iota(nto_ixs.begin(), nto_ixs.end(), 0);
        const size_t g_size = num_var * num_var;
        std::vector<int> G(g_size, 1);
        ReducedGC gc = {
            num_var,
            num_phen,
            args.max_level_one,
            nto_ixs, // new_to_old_indices
            G,
            sq_corrs,
            sq_ess
        };
        float th = hetcor_threshold(args.alpha);

        std::cout << "Starting first cusk stage" << std::endl;
        gc = run_cusk(
            gc,
            th,
            args.depth,
            args.max_level_one,
            time_index_traits
        );
        gc.to_file(make_path(outdir, "trait_only", ""));
    } else {
        std::cout << "Loading mxm" << std::endl;
        MarkerSummaryStats mxm = MarkerSummaryStats(mxm_path);
        
        std::cout << "Loading mxp summary stats" << std::endl;

        if (args.hetcor && args.merged) {
            MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, mxp_se_path, marker_ixs);
        } else if (args.hetcor) {
            MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, mxp_se_path, block);
        } else if (args.merged) {
            MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, marker_ixs);
        } else {
            MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, block);
        }

        // check if all dims check out
        if (pxp.get_num_phen() != mxp.get_num_phen())
        {
            std::cout << "Numbers of traits seem to differ between pxp and mxp" << std::endl;
            exit(1);
        }
        if (mxm.get_num_markers() != mxp.get_num_markers())
        {
            std::cout << "Numbers of markers seem to differ between mxm and mxp" << std::endl;
            std::cout << "mxp: " << mxp.get_num_markers() << " x " << mxp.get_num_phen() << std::endl;
            std::cout << "mxm: " << mxm.get_num_markers() << " x " << mxm.get_num_markers()
                    << std::endl;
            exit(1);
        }

        // dims ok, now merge them all into one matrix.
        std::cout << "Merging correlations into single matrix" << std::endl;
        size_t num_phen = pxp.get_num_phen();
        size_t num_markers = mxm.get_num_markers();
        size_t num_var = num_markers + num_phen;
        CuskssSquareInputs csi = make_square_mxm_mxp_pxp_corrmat(
            mxm,
            mxp,
            pxp,
            args.pearson_sample_size,
            args.hetcor // heterogoneous_sample_size ?
        );
        std::vector<float> sq_corrs = csi.correlations;
        std::vector<float> sq_ess = csi.sample_sizes;

        // init new-to-old indices for initial gc
        std::vector<int> nto_ixs(num_var);
        std::iota(v.begin(), v.end(), 0);
        const size_t g_size = num_var * num_var;
        std::vector<int> G(g_size, 1);
        ReducedGC gc = {
            num_var,
            num_phen,
            args.max_level_one,
            nto_ixs, // new_to_old_indices
            G,
            sq_corrs,
            sq_ess
        };
        float th = hetcor_threshold(alpha);

        std::cout << "Starting first cusk stage" << std::endl;
        gc = run_cusk(
            gc,
            th,
            args.depth,
            args.max_level_one,
            time_index_traits
        );

        if (args.two_stage) {
            std::cout << "Starting second cusk stage" << std::endl;
            gc = run_cusk(
                gc,
                th,
                args.depth,
                args.max_level_two,
                time_index_traits
            )
        }
        std::cout << "Retained " << gc.num_markers() << " markers" << std::endl;
        if (args.merged) {
            gc.to_file(make_path(outdir, "cuskss_merged", ""))
        } else {
            gc.to_file(make_path(outdir, block.to_file_string(), ""))
        }
    }
}

const std::string BLOCK_USAGE = R"(
Group markers into approximately unlinked blocks

usage: mps block <bfiles> <max-block-size> <device-mem-gb> <corr-width>

arguments:
    bfiles          stem of .bed, .bim, .fam files
    max-block-size  maximum number of markers per block
    device-mem-gb   maximum memory available on gpu in GB
    corr-width      max distance at which to compute correlations
)";

const int BLOCK_NARGS = 6;

void make_blocks(int argc, char *argv[])
{
    check_nargs(argc, BLOCK_NARGS, BLOCK_USAGE);

    std::string bed_base_path = argv[2];
    int max_block_size = std::stoi(argv[3]);
    float device_mem_gb = std::stof(argv[4]);
    size_t corr_width = std::stoi(argv[5]);

    std::cout << "Checking paths" << std::endl;
    check_bed_path(bed_base_path);

    BfilesBase bfiles(bed_base_path);
    BedDims dim(bfiles);
    BimInfo bim(bfiles.bim());

    size_t bytes_per_bed_col = ((dim.get_num_samples() + 3) / 4);

    for (auto cid : bim.chr_ids)
    {
        size_t num_markers = bim.get_num_markers_on_chr(cid);
        size_t mem_host_gb =
            ((num_markers * bytes_per_bed_col) + corr_width * num_markers * 4) * std::pow(10, -9);
        std::cout << "[Chr " << cid << "]: At least " << mem_host_gb
                  << " GB in host memory required." << std::endl;

        std::cout << "[Chr " << cid << "]: Loading bed data for " << num_markers << " markers."
                  << std::endl;
        std::vector<unsigned char> chr_bed = read_chr_from_bed(bfiles.bed(), cid, bim, dim);

        std::cout << "[Chr " << cid << "]: Computing correlations." << std::endl;
        std::vector<float> mcorrs = cal_mcorrk_banded(
            chr_bed, dim, bim.get_num_markers_on_chr(cid), corr_width, device_mem_gb
        );

        std::cout << "[Chr " << cid << "]: Computing row sums." << std::endl;
        std::vector<float> mcorrkb_row_sums =
            marker_corr_banded_mat_row_abs_sums(mcorrs, corr_width);

        std::cout << "[Chr " << cid << "]: Making blocks." << std::endl;
        std::vector<MarkerBlock> blocks = block_chr(mcorrkb_row_sums, cid, max_block_size);

        std::cout << "[Chr " << cid << "]: Partitioned into " << blocks.size() << " blocks."
                  << std::endl;
        std::cout << "[Chr " << cid << "]: Writing blocks to output file." << std::endl;
        write_marker_blocks_to_file(blocks, bfiles.blocks(max_block_size));
    }

    std::cout << "Done." << std::endl;
}

const std::string CUSK_SINGLE_USAGE = R"(
Run cuPC on a single block of a block diagonal genomic covariance matrix.

usage: mps cusk <.phen> <bfiles> <.blocks> <alpha> <max-level> <depth> <outdir> <first-block>

arguments:
    .phen                     path to standardized phenotype tsv
    bfiles                    stem of .bed, .means, .stds, .dim files
    .blocks                   file with genomic block definitions
    alpha                     significance level
    max-level-stage-one       maximal size of seperation sets in cuPC round one ( <= 14)
    max-level-stage_two       maximal size of seperation sets in cuPC round two ( <= 14)
    depth                     max depth at which marker variables are kept as ancestors
    outdir                    outdir
    block-index               0-based index of the block to run cupc for
)";

const int CUSK_SINGLE_NARGS = 11;

void cusk(int argc, char *argv[])
{
    check_nargs(argc, CUSK_SINGLE_NARGS, CUSK_SINGLE_USAGE);

    std::cout << "Got args: " << std::endl;
    std::string phen_path = argv[2];
    std::cout << ".phen: " << phen_path << std::endl;
    std::string bed_base_path = argv[3];
    std::cout << "bfiles: " << bed_base_path << std::endl;
    std::string block_path = argv[4];
    std::cout << ".blocks: " << block_path << std::endl;
    float alpha = std::stof(argv[5]);
    std::cout << "alpha: " << alpha << std::endl;
    int max_level = std::stoi(argv[6]);
    std::cout << "max_level: " << max_level << std::endl;
    int max_level_two = std::stoi(argv[7]);
    std::cout << "max_level: " << max_level << std::endl;
    int depth = std::stoi(argv[8]);
    std::cout << "depth: " << depth << std::endl;
    std::string outdir = (std::string)argv[9];
    std::cout << "outdir: " << outdir << std::endl;
    int block_index = std::stoi(argv[10]);
    std::cout << "block-index: " << block_index << std::endl;

    std::cout << "Checking paths" << std::endl;

    check_prepped_bed_path(bed_base_path);
    check_path(phen_path);
    check_path(block_path);
    check_path(outdir);

    Phen phen = load_phen(phen_path);
    BfilesBase bfiles(bed_base_path);
    BedDims dims(bfiles.dim());

    if (phen.get_num_samples() != dims.get_num_samples())
    {
        std::cerr << "different num samples in phen and dims" << std::endl;
        exit(1);
    }

    BimInfo bim(bfiles.bim());
    size_t num_individuals = dims.get_num_samples();
    size_t num_phen = phen.get_num_phen();
    size_t phen_corr_mat_size = num_phen * (num_phen - 1) / 2;
    std::vector<float> phen_corr(phen_corr_mat_size, 0.0);
    std::cout << "Found " << num_phen << " phenotypes" << std::endl;

    std::cout << "Loading blocks" << std::endl;
    std::vector<MarkerBlock> blocks = read_blocks_from_file(block_path);

    std::cout << "Found " << blocks.size() << " blocks" << std::endl;

    for (auto b : blocks)
    {
        if (b.get_first_marker_ix() >= bim.get_num_markers_on_chr(b.get_chr_id()) ||
            (b.get_last_marker_ix() >= bim.get_num_markers_on_chr(b.get_chr_id())))
        {
            std::cerr << "block out of bounds with first_ix: " << b.get_first_marker_ix()
                      << " last_ix: " << b.get_last_marker_ix() << std::endl;
            exit(1);
        }
    }

    // size_t num_individuals = dims.get_num_samples();
    std::vector<float> Th = threshold_array(num_individuals, alpha);

    std::cout << "Number of levels: " << max_level << std::endl;

    std::cout << "Setting level thr for cuPC: " << std::endl;
    for (int i = 0; i <= max_level; ++i)
    {
        std::cout << "\t Level: " << i << " thr: " << Th[i] << std::endl;
    }

    size_t bid = block_index;

    std::cout << std::endl;
    std::cout << "Processing block " << bid + 1 << " / " << blocks.size() << std::endl;

    MarkerBlock block = blocks[bid];

    size_t num_markers = block.block_size();

    std::cout << "Block size: " << num_markers << std::endl;

    std::cout << "Loading bed data" << std::endl;

    // load block data
    std::vector<unsigned char> bedblock = read_block_from_bed(bfiles.bed(), block, dims, bim);

    size_t chr_start = bim.get_global_chr_start(block.get_chr_id());
    std::vector<float> means = read_floats_from_line_range(
        bfiles.means(),
        chr_start + block.get_first_marker_ix(),
        chr_start + block.get_last_marker_ix()
    );
    std::vector<float> stds = read_floats_from_line_range(
        bfiles.stds(),
        chr_start + block.get_first_marker_ix(),
        chr_start + block.get_last_marker_ix()
    );

    if ((means.size() != block.block_size()) || (stds.size() != block.block_size()))
    {
        std::cout << "block size and number of means or stds differ" << std::endl;
        exit(1);
    }

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    size_t marker_phen_corr_mat_size = num_markers * num_phen;
    std::vector<float> marker_phen_corr(marker_phen_corr_mat_size, 0.0);

    std::cout << "Checking for significant marker - phen correlations" << std::endl;

    cu_marker_phen_corr_pearson(
        bedblock.data(),
        phen.data.data(),
        num_markers,
        num_individuals,
        num_phen,
        means.data(),
        stds.data(),
        marker_phen_corr.data()
    );

    int num_sig_corrs = 0;
    for (float c : marker_phen_corr)
    {
        num_sig_corrs += (fabs(0.5 * (log(fabs((1 + c))) - log(fabs(1 - c)))) >= Th[0]);
    }

    if (num_sig_corrs > 0)
    {
        std::cout << "Found " << num_sig_corrs << " marker - phen correlations. Proceeding."
                  << std::endl;
    }
    else
    {
        std::cout << "No significant correlations found. Skipping block." << std::endl;
        return;
    }

    std::cout << "Computing all correlations" << std::endl;

    // compute correlations
    cu_corr_pearson_npn(
        bedblock.data(),
        phen.data.data(),
        num_markers,
        num_individuals,
        num_phen,
        means.data(),
        stds.data(),
        marker_corr.data(),
        marker_phen_corr.data(),
        phen_corr.data()
    );

    std::cout << "Reformating corrs to n2 format" << std::endl;

    // make n2 matrix to please cuPC
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs(num_var * num_var, 1.0);

    size_t sq_row_ix = 0;
    size_t sq_col_ix = 1;
    for (size_t i = 0; i < marker_corr_mat_size; ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)marker_corr[i];
        sq_corrs[num_var * sq_col_ix + sq_row_ix] = (float)marker_corr[i];
        if (sq_col_ix == num_markers - 1)
        {
            ++sq_row_ix;
            sq_col_ix = sq_row_ix + 1;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    sq_row_ix = 0;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < marker_phen_corr_mat_size; ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)marker_phen_corr[i];
        sq_corrs[num_var * sq_col_ix + sq_row_ix] = (float)marker_phen_corr[i];
        if (sq_col_ix == (num_var - 1))
        {
            sq_col_ix = num_markers;
            ++sq_row_ix;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    sq_row_ix = num_markers;
    sq_col_ix = num_markers + 1;
    for (size_t i = 0; i < phen_corr_mat_size; ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)phen_corr[i];
        sq_corrs[num_var * sq_col_ix + sq_row_ix] = (float)phen_corr[i];
        if (sq_col_ix == (num_var - 1))
        {
            ++sq_row_ix;
            sq_col_ix = sq_row_ix + 1;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    if (WRITE_FULL_CORRMATS)
    {
        write_floats_to_binary(
            sq_corrs.data(),
            sq_corrs.size(),
            make_path(outdir, block.to_file_string(), ".all_corrs")
        );
    }

    std::cout << "Running cuPC" << std::endl;

    // call cuPC
    int p = num_var;
    const size_t sepset_size = p * p * ML;
    const size_t g_size = num_var * num_var;
    std::vector<float> pmax(g_size, 0.0);
    std::vector<int> G(g_size, 1);
    std::vector<int> sepset(sepset_size, 0);
    int l = 0;
    Skeleton(sq_corrs.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data());

    std::unordered_set<int> variable_subset = subset_variables(G, num_var, num_markers, depth);
    ReducedGCS gcs = reduce_gcs(G, sq_corrs, sepset, variable_subset, num_var, num_phen, max_level);
    std::cout << "Starting second cusk stage" << std::endl;
    gcs = reduced_gcs_cusk(gcs, Th, depth, max_level_two);
    std::cout << "Retained " << gcs.num_markers() << " markers" << std::endl;
    gcs.to_file(make_path(outdir, block.to_file_string(), ""));
}


const std::string PREP_USAGE = R"(
Prepare input (PLINK) .bed file for mps.

Compute column means and standard deviations, create all files necessary for mps bdpc.

usage: mps prep <.bfiles>

arguments:
    .bfiles filestem of .bed, .bim, .fam fileset
)";

const int PREP_NARGS = 3;

// TODO: this should also parse .phen files and make sure that no info is missing and order
// checks out
void prep_bed(int argc, char *argv[])
{
    // check if enough args present
    if ((argc != PREP_NARGS) || (argv[2] == "--help") || (argv[2] == "-h"))
    {
        std::cout << PREP_USAGE << std::endl;
        exit(1);
    }

    std::string bed_path = (std::string)argv[2];
    check_bed_path(bed_path);
    prep_bed_no_impute(BfilesBase(bed_path));
}


const std::string MPS_USAGE = R"(
usage: mps <command> [<args>]

commands:
    prep                    Prepare input (PLINK) .bed file for mps
    block                   Tile marker x marker correlation matrix
    cusk                    Run cuda-skeleton on a single block of block diagonal genomic covariance matrix
    cuskss                  Run cuda-skeleton on a block of markers and traits with pre-computed correlations.

contact:
    nick.machnik@gmail.com
)";

auto main(int argc, char *argv[]) -> int
{
    if (argc == 1)
    {
        std::cout << MPS_USAGE << std::endl;
        return EXIT_SUCCESS;
    }

    std::string cmd = (std::string)argv[1];

    if ((cmd == "--help") || (cmd == "-h"))
    {
        std::cout << MPS_USAGE << std::endl;
    }
    else if (cmd == "cuskss")
    {
        std::string mxm_path = argv[2];
        std::string mxp_path = argv[3];
        std::string mxp_se_path = argv[4];
        std::string pxp_path = argv[5];
        std::string pxp_se_path = argv[6];
        std::string time_index_path = argv[7];
        int block_index = std::stoi(argv[8]);
        std::string blockfile_path = argv[9];
        std::string marker_index_path = argv[10];
        float alpha = std::stof(argv[11]);
        int max_level_one = std::stoi(argv[12]);
        int max_level_two = std::stoi(argv[13]);
        int max_depth = std::stoi(argv[14]);
        int num_samples = std::stoi(argv[15]);
        std::string outdir_path = argv[16];

        bool merged = marker_index_path != "NULL";
        bool hetcor = mxp_se_path != "NULL";
        bool trait_only = mxm_path != "NULL";
        bool two_stage = max_level_two > 0;
        bool time_indexed = time_index_path != "NULL";

        // required args
        check_path(pxp_path);
        check_path(block_path);
        check_path(outdir);

        if (hetcor || merged || !trait_only) {
            check_path(mxm_path);
            check_path(mxp_path);
        }

        if (hetcor) {
            check_path(mxp_se_path);
            check_path(pxp_se_path);
        }

        if (time_indexed) {
            check_path(time_index_path);
        }

        CuskssArgs args = {
            merged,             // merged
            hetcor,             // hetcor
            trait_only,         // trait_only
            two_stage,          // two_stage
            time_indexed,       // time_indexed
            alpha,              // alpha
            num_samples,        // pearson_sample_size
            max_level_one,      // max_level
            max_level_two,      // max_level_two
            max_depth,          // depth
            block_index,        // block_ix
            blockfile_path,         // block_path
            marker_index_path,  // marker_ixs_path
            mxm_path,           // mxm_path
            mxp_path,           // mxp_path
            mxp_se_path,        // mxp_se_path
            pxp_path,           // pxp_path
            pxp_se_path,        // pxp_se_path
            time_index_path,    // time_index_path
            outdir_path         // outdir
        };
        cuskss(args);
    }
    else if (cmd == "prep")
    {
        prep_bed(argc, argv);
    }
    else if (cmd == "cusk")
    {
        cusk(argc, argv);
    }
    else if (cmd == "block")
    {
        make_blocks(argc, argv);
    }
    else
    {
        std::cout << MPS_USAGE << std::endl;
    }

    return EXIT_SUCCESS;
}
