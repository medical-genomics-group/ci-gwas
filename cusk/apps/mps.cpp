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

#include <array>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

bool WRITE_FULL_CORRMATS = false;

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

ReducedGC reduced_gc_cusk(
    ReducedGC gc, std::vector<float> &sample_sizes, float threshold, int max_depth, int max_level
)
{
    int num_var = gc.num_var;
    int start_level = 0;
    hetcor_skeleton(
        gc.C.data(),
        &num_var,
        gc.G.data(),
        sample_sizes.data(),
        &threshold,
        &start_level,
        &max_level
    );
    std::unordered_set<int> variable_subset =
        subset_variables(gc.G, gc.num_var, gc.num_markers(), max_depth);
    return reduce_gc(
        gc.G, gc.C, variable_subset, gc.num_var, gc.num_phen, ML, gc.new_to_old_indices
    );
}

ReducedGCS reduced_gcs_cusk(
    ReducedGCS gcs, std::vector<float> &thresholds, int max_depth, int max_level
)
{
    int num_var = gcs.num_var;
    int start_level = 0;
    const size_t sepset_size = gcs.num_var * gcs.num_var * ML;
    std::vector<int> sepset(sepset_size, 0);
    const size_t g_size = gcs.num_var * gcs.num_var;
    std::vector<float> pmax(g_size, 0.0);
    Skeleton(
        gcs.C.data(),
        &num_var,
        gcs.G.data(),
        thresholds.data(),
        &start_level,
        &max_level,
        pmax.data(),
        sepset.data()
    );
    std::unordered_set<int> variable_subset =
        subset_variables(gcs.G, gcs.num_var, gcs.num_markers(), max_depth);
    return reduce_gcs(
        gcs.G, gcs.C, sepset, variable_subset, gcs.num_var, gcs.num_phen, ML, gcs.new_to_old_indices
    );
}

const std::string MCORRKB_CHR_USAGE = R"(
Compute the banded Kendall correlation matrix for a given chromosome

usage: mps mcorrkb-chr <bfiles> <device-mem-gb> <corr-width>

arguments:
    bfiles          stem of .bed, .bim, .fam files
    chr             id of chr (as in .bim)
    device-mem-gb   maximum memory available on gpu in GB
    corr-width      max distance at which to compute correlations    
)";

const int MCORRKB_CHR_NARGS = 6;

void mcorrkb_chr(int argc, char *argv[])
{
    check_nargs(argc, MCORRKB_CHR_NARGS, MCORRKB_CHR_USAGE);

    std::string bed_base_path = argv[2];
    std::string chr_id = argv[3];
    float device_mem_gb = std::stof(argv[4]);
    size_t corr_width = std::stoi(argv[5]);

    std::cout << "Checking paths" << std::endl;
    check_bed_path(bed_base_path);

    BfilesBase bfiles(bed_base_path);
    BedDims dim(bfiles);
    BimInfo bim(bfiles.bim());

    size_t bytes_per_bed_col = ((dim.get_num_samples() + 3) / 4);

    size_t num_markers = bim.get_num_markers_on_chr(chr_id);
    size_t mem_host_gb =
        ((num_markers * bytes_per_bed_col) + corr_width * num_markers * 4) * std::pow(10, -9);
    std::cout << "[Chr " << chr_id << "]: At least " << mem_host_gb
              << " GB in host memory required." << std::endl;

    std::cout << "[Chr " << chr_id << "]: Loading bed data for " << num_markers << " markers."
              << std::endl;
    std::vector<unsigned char> chr_bed = read_chr_from_bed(bfiles.bed(), chr_id, bim, dim);

    std::cout << "[Chr " << chr_id << "]: Computing correlations." << std::endl;
    std::vector<float> mcorrs = cal_mcorrk_banded(
        chr_bed, dim, bim.get_num_markers_on_chr(chr_id), corr_width, device_mem_gb
    );

    write_floats_to_binary(mcorrs.data(), mcorrs.size(), bfiles.with_suffix(chr_id + ".mcorrkb"));
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

const std::string PHENOPC_USAGE = R"(
Run cuPC on phenotypes only.

usage: mps phenpc <.phen> <alpha> <max-level> <outdir>

arguments:
    .phen           path to standardized phenotype tsv
    alpha           significance level
    max-level       maximal size of seperation sets in cuPC ( <= 14)
    outdir          outdir
)";

const int PHENOPC_NARGS = 6;

void phenotype_pc(int argc, char *argv[])
{
    std::cout << "Computing skeleton for phenotypes only" << std::endl;
    check_nargs(argc, PHENOPC_NARGS, PHENOPC_USAGE);

    std::string phen_path = argv[2];
    float alpha = std::stof(argv[3]);
    int max_level = std::stoi(argv[4]);
    std::string outdir = (std::string)argv[5];

    Phen phen = load_phen(phen_path);
    size_t num_phen = phen.get_num_phen();
    size_t phen_corr_mat_size = num_phen * (num_phen - 1) / 2;
    std::vector<float> phen_corr(phen_corr_mat_size, 0.0);
    std::cout << "Found " << num_phen << " phenotypes" << std::endl;

    // size_t num_individuals = dims.get_num_samples();
    std::vector<float> Th = threshold_array(phen.get_num_samples(), alpha);

    std::cout << "Number of levels: " << max_level << std::endl;

    std::cout << "Setting level thr for cuPC: " << std::endl;
    for (int i = 0; i <= max_level; ++i)
    {
        std::cout << "\t Level: " << i << " thr: " << Th[i] << std::endl;
    }

    std::cout << "Computing correlations" << std::endl;
    cu_phen_corr_pearson_npn(
        phen.data.data(), phen.get_num_samples(), phen.get_num_phen(), phen_corr.data()
    );

    std::cout << "Reformating corrs to n2 format" << std::endl;
    std::vector<float> sq_corrs(num_phen * num_phen, 1.0);

    size_t sq_row_ix = 0;
    size_t sq_col_ix = 1;
    for (size_t i = 0; i < phen_corr_mat_size; ++i)
    {
        sq_corrs[num_phen * sq_row_ix + sq_col_ix] = phen_corr[i];
        sq_corrs[num_phen * sq_col_ix + sq_row_ix] = phen_corr[i];
        if (sq_col_ix == num_phen - 1)
        {
            ++sq_row_ix;
            sq_col_ix = sq_row_ix + 1;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    std::cout << "Running cuPC" << std::endl;

    // call cuPC
    int p = num_phen;
    const size_t sepset_size = p * p * ML;
    const size_t g_size = p * p;
    std::vector<float> pmax(g_size, 0.0);
    std::vector<int> G(g_size, 1);
    std::vector<int> sepset(sepset_size, 0);
    int l = 0;
    Skeleton(sq_corrs.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data());

    std::unordered_set<int> variable_subset = {};
    for (size_t i = 0; i < num_phen; i++)
    {
        variable_subset.insert(i);
    }

    ReducedGCS gcs = reduce_gcs(G, sq_corrs, sepset, variable_subset, p, num_phen, max_level);
    gcs.to_file(make_path(outdir, "pheno_sk", ""));
}

const std::string SIMPC_USAGE = R"(
Like bdpcss but with mtx file input and with a single block.

usage: mps simpc <mtx> <alpha> <max-level> <depth> <num_samples> <num-phen> <outdir>

arguments:
    mtx             correlation matrix where the last `num-phen` variables are the phenotypes. In .mtx format.
    alpha           significance level
    max-level       maximal size of seperation sets in cuPC
    depth           max depth at which marker variables are kept as ancestors
    num-samples     number of samples used in correlation calculations
    num-phen        number of phenotypes
    outdir          outdir
)";

const int SIMPC_NARGS = 9;

void sim_pc(int argc, char *argv[])
{
    check_nargs(argc, SIMPC_NARGS, SIMPC_USAGE);

    std::string mtx_path = argv[2];
    float alpha = std::stof(argv[3]);
    int max_level = std::stoi(argv[4]);
    int depth = std::stoi(argv[5]);
    int num_individuals = std::stoi(argv[6]);
    int num_phen = std::stoi(argv[7]);
    std::string outdir = (std::string)argv[8];

    check_path(mtx_path);

    std::cout << "Loading phenotype correlations" << std::endl;
    std::vector<float> corrs = read_correlations_from_mtx(mtx_path);

    int num_var = (int)floor(sqrt(corrs.size()));
    int num_markers = num_var - num_phen;

    // size_t num_individuals = dims.get_num_samples();
    std::vector<float> Th = threshold_array(num_individuals, alpha);

    std::cout << "Number of levels: " << max_level << std::endl;

    std::cout << "Setting level thr for cuPC: " << std::endl;
    for (int i = 0; i <= max_level; ++i)
    {
        std::cout << "\t Level: " << i << " thr: " << Th[i] << std::endl;
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
    Skeleton(corrs.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data());

    std::cout << "Reducing data to phenotype parent sets" << std::endl;

    std::unordered_set<int> variable_subset = subset_variables(G, num_var, num_markers, depth);

    ReducedGCS gcs = reduce_gcs(G, corrs, sepset, variable_subset, num_var, num_phen, max_level);

    std::cout << "Retained " << (variable_subset.size() - num_phen) << " / " << num_markers
              << " markers" << std::endl;

    gcs.to_file(make_path(outdir, "skeleton", ""));
}

const std::string MP_CORR_DIST_USAGE = R"(
Compute all marker-phenotype and phenotype-phenotype correlations.

usage: mps bdpc <.phen> <bfiles> <.blocks> <outdir>

arguments:
    .phen           path to standardized phenotype tsv
    bfiles          stem of .bed, .means, .stds, .dim files
    .blocks         file with genomic block definitions
    outdir          outdir
)";

const int MP_CORR_DIST_NARGS = 6;

void marker_pheno_corr_dist(int argc, char *argv[])
{
    check_nargs(argc, MP_CORR_DIST_NARGS, MP_CORR_DIST_USAGE);

    std::string phen_path = argv[2];
    std::string bed_base_path = argv[3];
    std::string block_path = argv[4];
    std::string outdir = (std::string)argv[5];

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

    std::cout << "Computing phenotype correlations" << std::endl;
    cu_phen_corr_pearson_npn(
        phen.data.data(), phen.get_num_samples(), phen.get_num_phen(), phen_corr.data()
    );

    std::cout << "Reformating corrs to n2 format" << std::endl;
    std::vector<float> sq_corrs(num_phen * num_phen, 1.0);

    size_t sq_row_ix = 0;
    size_t sq_col_ix = 1;
    for (size_t i = 0; i < phen_corr_mat_size; ++i)
    {
        sq_corrs[num_phen * sq_row_ix + sq_col_ix] = phen_corr[i];
        sq_corrs[num_phen * sq_col_ix + sq_row_ix] = phen_corr[i];
        if (sq_col_ix == num_phen - 1)
        {
            ++sq_row_ix;
            sq_col_ix = sq_row_ix + 1;
        }
        else
        {
            ++sq_col_ix;
        }
    }

    write_floats_to_binary(sq_corrs.data(), num_phen * num_phen, make_path(outdir, "pxp", ".corr"));

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

    for (size_t bid = 0; bid < blocks.size(); bid++)
    {
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

        size_t marker_phen_corr_mat_size = num_markers * num_phen;
        std::vector<float> marker_phen_corr(marker_phen_corr_mat_size, 0.0);

        std::cout << "Computing marker - phen correlations" << std::endl;

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

        append_floats_to_binary(
            marker_phen_corr.data(), marker_phen_corr_mat_size, make_path(outdir, "mxp", ".corr")
        );
    }
}

const std::string CUSKSS_TRAIT_ONLY_USAGE = R"(
Run cuda-skeleton on precomputed trait-trait correlations.

usage: mps cuskss-trait-only <pxp> <alpha> <max-level> <num-samples> <outdir>

arguments:
    pxp             correlations between all traits. Text format, rectangular, only upper triangle is used.
    alpha           significance level
    max-level       maximal size of seperation sets in cuPC ( <= 14)
    num-samples     number of samples used for computing correlations
    outdir          outdir
)";

const int CUSKSS_TRAIT_ONLY = 7;

void cuskss_trait_only(int argc, char *argv[])
{
    check_nargs(argc, CUSKSS_TRAIT_ONLY, CUSKSS_TRAIT_ONLY_USAGE);

    std::string pxp_path = argv[2];
    float alpha = std::stof(argv[3]);
    int max_level = std::stoi(argv[4]);
    int num_individuals = std::stoi(argv[5]);
    std::string outdir = (std::string)argv[6];

    check_path(pxp_path);
    check_path(outdir);

    // load everything
    std::cout << "Loading input files" << std::endl;
    std::cout << "Loading pxp" << std::endl;
    TraitSummaryStats pxp = TraitSummaryStats(pxp_path);

    // make n2 matrix to please cuPC
    size_t num_phen = pxp.get_num_phen();
    size_t num_markers = 0;
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs = pxp.get_corrs();

    std::vector<float> Th = threshold_array(num_individuals, alpha);

    // call cuPC
    int p = num_var;
    const size_t sepset_size = p * p * ML;
    const size_t g_size = num_var * num_var;
    std::vector<float> pmax(g_size, 0.0);
    std::vector<int> G(g_size, 1);
    std::vector<int> sepset(sepset_size, 0);
    int l = 0;
    Skeleton(sq_corrs.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data());

    std::unordered_set<int> variable_subset = {};
    for (size_t i = 0; i < num_phen; i++)
    {
        variable_subset.insert(i);
    }

    ReducedGCS gcs = reduce_gcs(G, sq_corrs, sepset, variable_subset, p, num_phen, max_level);
    gcs.to_file(make_path(outdir, "pheno_sk", ""));
}

const std::string CUSKSS_MERGED_USAGE = R"(
Run cuda-skeleton on pre-selected markers and traits with pre-computed correlations.

usage: mps cuskss-merged <mxm> <mxp> <pxp> <marker-ixs> <alpha> <max-level-stage-one> <max-level-stage-two> <depth> <num-samples> <outdir>

arguments:
    mxm                     correlations between markers in block. Binary of floats, lower triangular, with diagonal, row major.
    mxp                     correlations between markers in all blocks and all traits. Text format, rectangular.
    pxp                     correlations between all traits. Text format, rectangular, only upper triangle is used.
    marker-ixs              binary file of ints giving the row indices of the selected markers in the mxp file.
    alpha                   significance level
    max-level-stage-one     maximal size of seperation sets in cuPC round one ( <= 14)
    max-level-stage-two     maximal size of seperation sets in cuPC round two ( <= 14)
    depth                   max depth at which marker variables are kept as ancestors
    num-samples             number of samples used for computing correlations
    outdir                  outdir
)";

const int CUSKSS_MERGED_NARGS = 12;

void cuda_skeleton_summary_stats_merged_blocks(int argc, char *argv[])
{
    check_nargs(argc, CUSKSS_MERGED_NARGS, CUSKSS_MERGED_USAGE);

    std::string mxm_path = argv[2];
    std::string mxp_path = argv[3];
    std::string pxp_path = argv[4];
    std::string marker_ixs_path = argv[5];
    float alpha = std::stof(argv[6]);
    int max_level = std::stoi(argv[7]);
    int max_level_two = std::stoi(argv[8]);
    int depth = std::stoi(argv[9]);
    int num_individuals = std::stoi(argv[10]);
    std::string outdir = (std::string)argv[11];

    check_path(mxm_path);
    check_path(mxp_path);
    check_path(pxp_path);
    check_path(marker_ixs_path);
    check_path(outdir);

    // load everything
    std::cout << "Loading input files" << std::endl;
    std::cout << "Loading blocks" << std::endl;
    std::vector<int> marker_ixs = read_ints_from_binary(marker_ixs_path);
    std::cout << "Loading pxp" << std::endl;
    TraitSummaryStats pxp = TraitSummaryStats(pxp_path);
    std::cout << "Loading mxm" << std::endl;
    MarkerSummaryStats mxm = MarkerSummaryStats(mxm_path);
    std::cout << "Loading mxp" << std::endl;
    MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, marker_ixs);

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

    // make n2 matrix to please cuPC
    size_t num_phen = pxp.get_num_phen();
    size_t num_markers = mxm.get_num_markers();
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs(num_var * num_var, 1.0);
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
    sq_row_ix = 0;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < marker_phen_corr.size(); ++i)
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

    std::vector<float> phen_corr = pxp.get_corrs();
    sq_row_ix = num_markers;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < phen_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)phen_corr[i];
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

    if (WRITE_FULL_CORRMATS)
    {
        write_floats_to_binary(
            sq_corrs.data(), sq_corrs.size(), make_path(outdir, "cuskss_merged", ".all_corrs")
        );
    }

    std::cout << "Number of levels: " << max_level << std::endl;

    std::vector<float> Th = threshold_array(num_individuals, alpha);

    std::cout << "Setting level thr for cuPC: " << std::endl;
    for (int i = 0; i <= max_level; ++i)
    {
        std::cout << "\t Level: " << i << " thr: " << Th[i] << std::endl;
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
    gcs.to_file(make_path(outdir, "cuskss_merged", ""));
}

const std::string CUSKSS_USAGE = R"(
Run cuda-skeleton on a block of markers and traits with pre-computed correlations.

usage: mps cuskss <mxm> <mxp> <pxp> <block-ix> <.blocks> <alpha> <max-level> <depth> <num-samples> <outdir>

arguments:
    mxm                     correlations between markers in block. Binary of floats, lower triangular, with diagonal, row major.
    mxp                     correlations between markers in all blocks and all traits. Text format, rectangular.
    pxp                     correlations between all traits. Text format, rectangular, only upper triangle is used.
    block-ix                0-based index of the marker block in the .blocks file
    .blocks                 file with genomic block definitions
    alpha                   significance level
    max-level-stage-one     maximal size of seperation sets in cuPC round one ( <= 14)
    max-level-stage-two     maximal size of seperation sets in cuPC round two ( <= 14)
    depth                   max depth at which marker variables are kept as ancestors
    num-samples             number of samples used for computing correlations
    outdir                  outdir
)";

const int CUSKSS_NARGS = 13;

void cuda_skeleton_summary_stats(int argc, char *argv[])
{
    check_nargs(argc, CUSKSS_NARGS, CUSKSS_USAGE);

    std::string mxm_path = argv[2];
    std::string mxp_path = argv[3];
    std::string pxp_path = argv[4];
    int block_ix = std::stoi(argv[5]);
    std::string block_path = argv[6];
    float alpha = std::stof(argv[7]);
    int max_level = std::stoi(argv[8]);
    int max_level_two = std::stoi(argv[9]);
    int depth = std::stoi(argv[10]);
    int num_individuals = std::stoi(argv[11]);
    std::string outdir = (std::string)argv[12];

    check_path(mxm_path);
    check_path(mxp_path);
    check_path(pxp_path);
    check_path(block_path);
    check_path(outdir);

    // load everything
    std::cout << "Loading input files" << std::endl;
    std::cout << "Loading blocks" << std::endl;
    std::vector<MarkerBlock> blocks = read_blocks_from_file(block_path);
    MarkerBlock block = blocks[block_ix];
    std::cout << "Loading pxp" << std::endl;
    TraitSummaryStats pxp = TraitSummaryStats(pxp_path);
    std::cout << "Loading mxm" << std::endl;
    MarkerSummaryStats mxm = MarkerSummaryStats(mxm_path);
    std::cout << "Loading mxp" << std::endl;
    MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, block);

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

    // make n2 matrix to please cuPC
    size_t num_phen = pxp.get_num_phen();
    size_t num_markers = mxm.get_num_markers();
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs(num_var * num_var, 1.0);
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
    sq_row_ix = 0;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < marker_phen_corr.size(); ++i)
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

    std::vector<float> phen_corr = pxp.get_corrs();
    sq_row_ix = num_markers;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < phen_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)phen_corr[i];
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

    if (WRITE_FULL_CORRMATS)
    {
        write_floats_to_binary(
            sq_corrs.data(),
            sq_corrs.size(),
            make_path(outdir, block.to_file_string(), ".all_corrs")
        );
    }

    std::cout << "Number of levels: " << max_level << std::endl;

    std::vector<float> Th = threshold_array(num_individuals, alpha);

    std::cout << "Setting level thr for cuPC: " << std::endl;
    for (int i = 0; i <= max_level; ++i)
    {
        std::cout << "\t Level: " << i << " thr: " << Th[i] << std::endl;
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


const std::string CUSKSS_HET_USAGE = R"(
Run cuda-skeleton on a block of markers and traits with pre-computed heterogeneous correlations.

usage: mps cuskss <mxm> <mxp> <pxp> <block-ix> <.blocks> <alpha> <max-level> <depth> <num-samples> <outdir>

arguments:
    mxm                     correlations between markers in block. Binary of floats, lower triangular, with diagonal, row major.
    mxp                     correlations between markers in all blocks and all traits. Text format, rectangular.
    pxp                     correlations between all traits. Text format, rectangular, only upper triangle is used.
    block-ix                0-based index of the marker block in the .blocks file
    .blocks                 file with genomic block definitions
    alpha                   significance level
    max-level-stage-one     maximal size of seperation sets in cuPC round one ( <= 14)
    max-level-stage-two     maximal size of seperation sets in cuPC round two ( <= 14)
    depth                   max depth at which marker variables are kept as ancestors
    num-samples             number of samples used for computing correlations (square matrix of floats, num var * num_var entries, binary)
    outdir                  outdir
)";

const int CUSKSS_HET_NARGS = 13;

void cuda_skeleton_summary_stats_hetcor(int argc, char *argv[])
{
    check_nargs(argc, CUSKSS_HET_NARGS, CUSKSS_HET_USAGE);

    std::string mxm_path = argv[2];
    std::string mxp_path = argv[3];
    std::string pxp_path = argv[4];
    int block_ix = std::stoi(argv[5]);
    std::string block_path = argv[6];
    float alpha = std::stof(argv[7]);
    int max_level = std::stoi(argv[8]);
    int max_level_two = std::stoi(argv[9]);
    int depth = std::stoi(argv[10]);
    std::string ess_path = argv[11];
    std::string outdir = (std::string)argv[12];

    check_path(mxm_path);
    check_path(mxp_path);
    check_path(pxp_path);
    check_path(ess_path);
    check_path(block_path);
    check_path(outdir);

    // load everything
    std::cout << "Loading input files" << std::endl;
    std::cout << "Loading blocks" << std::endl;
    std::vector<MarkerBlock> blocks = read_blocks_from_file(block_path);
    MarkerBlock block = blocks[block_ix];
    std::cout << "Loading pxp" << std::endl;
    TraitSummaryStats pxp = TraitSummaryStats(pxp_path);
    std::cout << "Loading mxm" << std::endl;
    MarkerSummaryStats mxm = MarkerSummaryStats(mxm_path);
    std::cout << "Loading mxp" << std::endl;
    MarkerTraitSummaryStats mxp = MarkerTraitSummaryStats(mxp_path, block);
    std::cout << "Loading ess matrix" << std::endl;
    std::vector<float> ess = read_floats_from_binary(ess_path);

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
    if (ess.size() != (pxp.get_num_phen() + mxm.get_num_markers())) {
        std::cout << "Numbers of variables seem to differ between sample size matrix and pxp + mxm" << std::endl;
        exit(1);
    }

    // dims ok, now merge them all into one matrix.

    std::cout << "Merging correlations into single matrix" << std::endl;

    // make n2 matrix to please cuPC
    size_t num_phen = pxp.get_num_phen();
    size_t num_markers = mxm.get_num_markers();
    size_t num_var = num_markers + num_phen;
    std::vector<float> sq_corrs(num_var * num_var, 1.0);
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
    sq_row_ix = 0;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < marker_phen_corr.size(); ++i)
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

    std::vector<float> phen_corr = pxp.get_corrs();
    sq_row_ix = num_markers;
    sq_col_ix = num_markers;
    for (size_t i = 0; i < phen_corr.size(); ++i)
    {
        sq_corrs[num_var * sq_row_ix + sq_col_ix] = (float)phen_corr[i];
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

    if (WRITE_FULL_CORRMATS)
    {
        write_floats_to_binary(
            sq_corrs.data(),
            sq_corrs.size(),
            make_path(outdir, block.to_file_string(), ".all_corrs")
        );
    }

    float th = hetcor_threshold(alpha);

    std::cout << "Running cuPC" << std::endl;

    // call cuPC
    int p = num_var;
    const size_t g_size = num_var * num_var;
    std::vector<int> G(g_size, 1);
    int init_level = 0;
    hetcor_skeleton(sq_corrs.data(), &p, G.data(), ess.data(), &th, &init_level, &max_level);

    std::unordered_set<int> variable_subset = subset_variables(G, num_var, num_markers, depth);
    ReducedGC gc = reduce_gc(G, sq_corrs, variable_subset, num_var, num_phen, max_level);
    std::cout << "Starting second cusk stage" << std::endl;
    gc = reduced_gc_cusk(gc, ess, th, depth, max_level_two);
    std::cout << "Retained " << gc.num_markers() << " markers" << std::endl;
    gc.to_file(make_path(outdir, block.to_file_string(), ""));
}

const std::string CUSK_USAGE = R"(
Run cuPC on block diagonal genomic covariance matrix.

usage: mps cusk <.phen> <bfiles> <.blocks> <alpha> <max-level> <depth> <outdir> <first-block>

arguments:
    .phen                       path to standardized phenotype tsv
    bfiles                      stem of .bed, .means, .stds, .dim files
    .blocks                     file with genomic block definitions
    alpha                       significance level
    max-level-stage-one         maximal size of seperation sets in cuPC round one ( <= 14)
    max-level-stage-two         maximal size of seperation sets in cuPC round two ( <= 14)
    depth                       max depth at which marker variables are kept as ancestors
    outdir                      outdir
    first-block                 0-based index first block to start cusk from. Meant for resuming interrupted jobs.
)";

const int CUSK_NARGS = 11;

void cusk(int argc, char *argv[])
{
    check_nargs(argc, CUSK_NARGS, CUSK_USAGE);

    std::string phen_path = argv[2];
    std::string bed_base_path = argv[3];
    std::string block_path = argv[4];
    float alpha = std::stof(argv[5]);
    int max_level = std::stoi(argv[6]);
    int max_level_two = std::stoi(argv[7]);
    int depth = std::stoi(argv[8]);
    std::string outdir = (std::string)argv[9];
    int first_block = std::stoi(argv[10]);

    std::cout << "Checking paths" << std::endl;

    check_prepped_bed_path(bed_base_path);
    check_path(phen_path);
    check_path(block_path);
    check_path(outdir);

    // call phenopc to compute skeleton phenotypes only
    char *phenopc_argv[PHENOPC_NARGS] = {argv[0], argv[1], argv[2], argv[5], argv[6], argv[8]};
    phenotype_pc(PHENOPC_NARGS, phenopc_argv);

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
    // TODO: testcase for num_phen = 1;
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

    for (size_t bid = first_block; bid < blocks.size(); bid++)
    {
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
            continue;
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

        std::cout << "Running cuPC" << std::endl;

        // call cuPC
        int p = num_var;
        const size_t sepset_size = p * p * ML;
        const size_t g_size = num_var * num_var;
        std::vector<float> pmax(g_size, 0.0);
        std::vector<int> G(g_size, 1);
        std::vector<int> sepset(sepset_size, 0);
        int l = 0;
        Skeleton(
            sq_corrs.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset.data()
        );

        std::unordered_set<int> variable_subset = subset_variables(G, num_var, num_markers, depth);
        ReducedGCS gcs =
            reduce_gcs(G, sq_corrs, sepset, variable_subset, num_var, num_phen, max_level);
        std::cout << "Starting second cusk stage" << std::endl;
        gcs = reduced_gcs_cusk(gcs, Th, depth, max_level_two);
        std::cout << "Retained " << gcs.num_markers() << " markers" << std::endl;
        gcs.to_file(make_path(outdir, block.to_file_string(), ""));
    }
}

const std::string CUSK_SINGLE_USAGE = R"(
Run cuPC on a single block of a block diagonal genomic covariance matrix.

usage: mps cusk-single <.phen> <bfiles> <.blocks> <alpha> <max-level> <depth> <outdir> <first-block>

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

void cusk_single(int argc, char *argv[])
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

const std::string ANTIDIAGSUMS_USAGE = R"(
Compute sums of anti-diagonals of marker correlation matrix.

usage: mps anti-diag-sums <corrs>

arguments:
    corrs    marker correlation matrix computes with `mps mcorrk` or `mps mcorrp`
)";

const int ANTIDIAGSUMS_NARGS = 3;

void anti_diag_sums(int argc, char *argv[])
{
    if (argc < ANTIDIAGSUMS_NARGS)
    {
        std::cout << ANTIDIAGSUMS_USAGE << std::endl;
        exit(1);
    }

    // check that path is valid
    std::string fpath = (std::string)argv[2];
    if (!path_exists(fpath))
    {
        std::cout << "file or directory not found: " << fpath << std::endl;
        exit(1);
    }

    printf("Loading correlations\n");
    fflush(stdout);
    std::vector<float> corrs = read_floats_from_binary(fpath);
    // get num markers from size of array

    size_t num_markers = num_variables_from_matrix_size(corrs.size());
    printf("Corr matrix has %i markers \n", num_markers);
    fflush(stdout);
    size_t num_anti_diagonals = 2 * num_markers - 3;

    printf("Allocating output vector\n");
    fflush(stdout);
    std::vector<float> sums(num_anti_diagonals, 0.0);

    printf("Computing antidiag sums\n");
    fflush(stdout);
    marker_corr_mat_antidiag_sums(num_markers, corrs.data(), sums.data());

    std::filesystem::path path(fpath);
    std::string dir = path.parent_path().string();
    std::string file = path.filename().string();

    write_single_column_file_with_suffix(sums, dir, file, ".ads");
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

const std::string MCORRK_USAGE = R"(
Compute Pearson correlation between markers as sin(pi / 2 * tau-b) where tau is the Kendall's tau-b

usage: mps mcorrk <prepdir> <chr> <device_mem_gb>

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dim files for each chromosome
    chr ID of the chromsome to be processed
    device_mem_gb   Amount of memory available on the GPU
)";

const int MCORRK_NARGS = 5;

void mcorrk(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < MCORRK_NARGS)
    {
        std::cout << MCORRK_USAGE << std::endl;
        exit(1);
    }

    float device_mem_gb = atof(argv[4]);
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[2] = {".bed", ".dim"};
    for (size_t i = 0; i < 2; ++i)
    {
        std::string fpath = make_path(out_dir, chr_id, req_suffixes[i]);
        if (!path_exists(fpath))
        {
            std::cout << "file or directory not found: " << fpath << std::endl;
            exit(1);
        }
    }

    // load data + check that file contents are valid
    // .dim
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dim");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dim file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    printf("Loading .bed \n");
    fflush(stdout);

    // .bed
    std::string bed_path = make_path(out_dir, chr_id, ".bed");
    size_t nbytes_per_block = (num_individuals + 3) / 4;
    size_t nbytes_marker_vals = num_markers * nbytes_per_block;
    std::vector<unsigned char> marker_vals(nbytes_marker_vals);
    read_n_bytes_from_binary(bed_path, nbytes_marker_vals, marker_vals);

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    // mem required for non-batched processing
    size_t nbytes_marker_data = num_individuals / 4 * num_markers;
    size_t nbytes_marker_corrs = num_markers * (num_markers - 1) * 2;
    size_t req_mem_bytes = nbytes_marker_data + nbytes_marker_corrs;

    if (req_mem_bytes > device_mem_bytes)
    {
        // figure out batch size
        float b = num_individuals / 4 - 2;
        size_t max_batch_size = (size_t)((-b + std::sqrt(b * b + 8 * device_mem_bytes)) / 4);
        size_t batch_nrows = max_batch_size / 2;

        printf("req mem: %zu, device mem: %zu \n", req_mem_bytes, device_mem_bytes);
        size_t mb_bytes =
            num_individuals / 4 * max_batch_size + max_batch_size * (max_batch_size - 1) * 2;
        printf("max_batch_size: %zu (%zu bytes \n)", max_batch_size, mb_bytes);

        printf("Device mem < required mem; Running tiled routine. \n");
        fflush(stdout);

        // printf("arg1: %p \n", marker_vals.data());
        // printf("arg2: %i \n", num_markers);
        // printf("arg3: %i \n", num_individuals);
        // printf("arg4: %i \n", batch_nrows);
        // printf("arg5: %p \n", marker_corr.data());
        // printf("wtf1 \n");
        // printf(
        //     "args: %p, %i, %i, %i, %p \n",
        //     marker_vals.data(),
        //     num_markers,
        //     num_individuals,
        //     batch_nrows,
        //     marker_corr.data()
        // );
        // printf("wtf2 \n");
        // printf("wtf3 \n");
        cu_marker_corr_pearson_npn_batched(
            marker_vals.data(), num_markers, num_individuals, batch_nrows, marker_corr.data()
        );
    }
    else
    {
        printf("Calling correlation main\n");
        fflush(stdout);
        // compute correlations
        cu_marker_corr_pearson_npn(
            marker_vals.data(), num_markers, num_individuals, marker_corr.data()
        );
    }

    printf("Writing results\n");
    fflush(stdout);
    // write results
    write_floats_to_binary(
        marker_corr.data(), marker_corr_mat_size, make_path(out_dir, chr_id, "_marker_corrk.bin")
    );

    std::cout << "Correlation matrix written to: " << out_dir << std::endl;
}

const std::string MCORRP_USAGE = R"(
Compute Pearson correlation between markers

usage: mps mcorrp <prepdir> <chr> <device_mem_gb>

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dim files for each chromosome
    chr ID of the chromsome to be processed
    device_mem_gb   Amount of memory available on the GPU
)";

const int MCORRP_NARGS = 5;

void mcorrp(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < MCORRP_NARGS)
    {
        std::cout << MCORRP_USAGE << std::endl;
        exit(1);
    }

    float device_mem_gb = atof(argv[4]);
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dim", ".means", ".stds"};
    for (size_t i = 0; i < 4; ++i)
    {
        std::string fpath = make_path(out_dir, chr_id, req_suffixes[i]);
        if (!path_exists(fpath))
        {
            std::cout << "file or directory not found: " << fpath << std::endl;
            exit(1);
        }
    }

    // load data + check that file contents are valid
    // .dim
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dim");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dim file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    printf("Loading .bed \n");
    fflush(stdout);

    // .bed
    std::string bed_path = make_path(out_dir, chr_id, ".bed");
    size_t nbytes_per_block = (num_individuals + 3) / 4;
    size_t nbytes_marker_vals = num_markers * nbytes_per_block;
    std::vector<unsigned char> marker_vals(nbytes_marker_vals);
    read_n_bytes_from_binary(bed_path, nbytes_marker_vals, marker_vals);

    printf("Loading means \n");
    fflush(stdout);

    // .means
    std::string means_path = make_path(out_dir, chr_id, ".means");
    std::vector<float> marker_means = read_floats_from_lines(means_path);
    assert((marker_means.size() == num_markers) && "number of marker means != num markers");

    printf("Loading stds \n");
    fflush(stdout);

    // .stds
    std::string stds_path = make_path(out_dir, chr_id, ".stds");
    std::vector<float> marker_stds = read_floats_from_lines(stds_path);
    assert((marker_stds.size() == num_markers) && "number of marker stds != num markers");

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    // mem required for non-batched processing
    size_t nbytes_marker_data = num_individuals / 4 * num_markers;
    size_t nbytes_marker_corrs = num_markers * (num_markers - 1) * 2;
    size_t req_mem_bytes = nbytes_marker_data + nbytes_marker_corrs;

    if (req_mem_bytes > device_mem_bytes)
    {
        // figure out batch size
        float b = num_individuals / 4 - 2;
        size_t max_batch_size = (size_t)((-b + std::sqrt(b * b + 8 * device_mem_bytes)) / 4);
        size_t batch_nrows = max_batch_size / 2;

        printf("Device mem < required mem; Running tiled routine. \n");
        fflush(stdout);

        cu_marker_corr_pearson_batched(
            marker_vals.data(),
            num_markers,
            num_individuals,
            marker_means.data(),
            marker_stds.data(),
            batch_nrows,
            marker_corr.data()
        );
    }
    else
    {
        printf("Calling correlation main\n");
        fflush(stdout);
        // compute correlations
        cu_marker_corr_pearson(
            marker_vals.data(),
            num_markers,
            num_individuals,
            marker_means.data(),
            marker_stds.data(),
            marker_corr.data()
        );
    }

    printf("Writing results\n");
    fflush(stdout);
    // write results
    write_floats_to_binary(
        marker_corr.data(), marker_corr_mat_size, make_path(out_dir, chr_id, "_marker_corrp.bin")
    );

    std::cout << "Correlation matrix written to: " << out_dir << std::endl;
}

const std::string SCORR_USAGE = R"(
Compute correlations between markers and phenotypes. The marker-marker correlation part 
is computed as a banded matrix.

usage: mps scorr <prepdir> <chr> <device_mem_gb> <distance_threshold> <.phen>

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dim files for each chromosome
    chr                 ID of the chromsome to be processed
    device_mem_gb       Amount of memory available on the GPU
    distance_threshold  Maximal distance between markers up to which correlations are computed (in bases)
    .phen               Path to .phen file with phenotype values, sorted in the same way as genotype info in the original bed file
)";

const int SCORR_NARGS = 7;

void scorr(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < SCORR_NARGS)
    {
        std::cout << SCORR_USAGE << std::endl;
        exit(1);
    }

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[5] = {".bed", ".dim", ".means", ".stds", ".bim"};
    for (size_t i = 0; i < 5; ++i)
    {
        std::string fpath = make_path(out_dir, chr_id, req_suffixes[i]);
        if (!path_exists(fpath))
        {
            std::cout << "file or directory not found: " << fpath << std::endl;
            exit(1);
        }
    }
    float device_mem_gb = atof(argv[4]);
    unsigned long dthr = atol(argv[5]);
    size_t corr_width = num_markers_within_distance(make_path(out_dir, chr_id, ".bim"), dthr);
    printf("Determined corr width: %u \n", corr_width);

    // .dim
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dim");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dim file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    printf("Loading .phen \n");
    // load phen data
    Phen phen = load_phen((std::string)argv[6]);
    size_t num_phen = phen.num_phen;

    assert(
        (phen.num_samples == num_individuals) &&
        "number of phen values != number of individuals in dim"
    );

    printf("num_individuals: %u \n", num_individuals);
    printf("num_phen: %u \n", num_phen);
    printf("num_markers: %u \n", num_markers);

    float device_mem_bytes = device_mem_gb * std::pow(10, 9);
    size_t upper_bound_batch_size = std::floor(
        (device_mem_bytes / 4 - (num_phen * num_individuals)) /
        (corr_width + num_phen + num_individuals / 16)
    );
    if (upper_bound_batch_size > num_markers)
    {
        upper_bound_batch_size = num_markers;
    }
    if (upper_bound_batch_size < corr_width)
    {
        printf(
            "Maximal batch size (%u) < corr width (%u). Decrease distance threshold or "
            "increase "
            "device memory. \n",
            upper_bound_batch_size,
            corr_width
        );
        exit(1);
    }

    printf("batch_size: %u \n", upper_bound_batch_size);

    printf("Loading .bed \n");

    // .bed
    std::string bed_path = make_path(out_dir, chr_id, ".bed");
    size_t nbytes_per_block = (num_individuals + 3) / 4;
    size_t nbytes_marker_vals = num_markers * nbytes_per_block;
    std::vector<unsigned char> marker_vals(nbytes_marker_vals);
    read_n_bytes_from_binary(bed_path, nbytes_marker_vals, marker_vals);

    printf("Loading means \n");

    // .means
    std::string means_path = make_path(out_dir, chr_id, ".means");
    std::vector<float> marker_means = read_floats_from_lines(means_path);
    assert((marker_means.size() == num_markers) && "number of marker means != num markers");

    printf("Loading stds \n");

    // .stds
    std::string stds_path = make_path(out_dir, chr_id, ".stds");
    std::vector<float> marker_stds = read_floats_from_lines(stds_path);
    assert((marker_stds.size() == num_markers) && "number of marker stds != num markers");

    printf("Allocating results array \n");

    // allocate correlation result arrays
    size_t corr_mat_size = (corr_width + num_phen) * (num_markers + num_phen);
    std::vector<float> corrs(corr_mat_size, 0.0);

    printf("Computing correlations\n");
    fflush(stdout);
    // compute correlations
    cu_corr_pearson_npn_batched_sparse(
        marker_vals.data(),
        phen.data.data(),
        num_markers,
        num_individuals,
        num_phen,
        marker_means.data(),
        marker_stds.data(),
        corr_width,
        upper_bound_batch_size,
        corrs.data()
    );

    printf("Writing results\n");
    // write results
    write_floats_to_binary(corrs.data(), corr_mat_size, make_path(out_dir, chr_id, "_scorr.bin"));

    std::cout << "Correlation matrix written to: " << out_dir << std::endl;
}

const std::string CORR_USAGE = R"(
Compute correlations between markers and phenotypes.

usage: mps corr <prepdir> <chr> <device_mem_gb> <.phen>...

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dim files for each chromosome
    chr ID of the chromsome to be processed
    device_mem_gb   Amount of memory available on the GPU
    .phen   Path to .phen file with phenotype values, sorted in the same way as genotype info in the original bed file
)";

const int CORR_NARGS = 6;

void corr(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < CORR_NARGS)
    {
        std::cout << CORR_USAGE << std::endl;
        exit(1);
    }

    float device_mem_gb = atof(argv[4]);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dim", ".means", ".stds"};
    for (size_t i = 0; i < 4; ++i)
    {
        std::string fpath = make_path(out_dir, chr_id, req_suffixes[i]);
        if (!path_exists(fpath))
        {
            std::cout << "file or directory not found: " << fpath << std::endl;
            exit(1);
        }
    }

    std::vector<std::string> phen_paths = {};
    for (size_t i = 5; i < argc; ++i)
    {
        std::string phen_path = (std::string)argv[i];
        if (!path_exists(phen_path))
        {
            std::cout << "file or directory not found: " << phen_path << std::endl;
            exit(1);
        }
        phen_paths.push_back(phen_path);
    }

    // load data + check that file contents are valid
    // .dim
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dim");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dim file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    // .phen
    size_t num_phen = phen_paths.size();
    std::vector<float> phen_vals = {};
    for (size_t i = 0; i < num_phen; ++i)
    {
        read_floats_from_lines(phen_paths[i], phen_vals);
    }

    printf("Loading .bed \n");

    // .bed
    std::string bed_path = make_path(out_dir, chr_id, ".bed");
    size_t nbytes_per_block = (num_individuals + 3) / 4;
    size_t nbytes_marker_vals = num_markers * nbytes_per_block;
    std::vector<unsigned char> marker_vals(nbytes_marker_vals);
    read_n_bytes_from_binary(bed_path, nbytes_marker_vals, marker_vals);

    printf("Loading means \n");

    // .means
    std::string means_path = make_path(out_dir, chr_id, ".means");
    std::vector<float> marker_means = read_floats_from_lines(means_path);
    assert((marker_means.size() == num_markers) && "number of marker means != num markers");

    printf("Loading stds \n");

    // .stds
    std::string stds_path = make_path(out_dir, chr_id, ".stds");
    std::vector<float> marker_stds = read_floats_from_lines(stds_path);
    assert((marker_stds.size() == num_markers) && "number of marker stds != num markers");

    printf("Allocating results arrays \n");

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    size_t marker_phen_corr_mat_size = num_markers * num_phen;
    std::vector<float> marker_phen_corr(marker_phen_corr_mat_size, 0.0);

    // TODO: testcase for num_phen = 1;
    size_t phen_corr_mat_size = num_phen * (num_phen - 1) / 2;
    std::vector<float> phen_corr(phen_corr_mat_size, 0.0);

    // mem required for non-batched processing
    size_t nbytes_marker_data = num_individuals / 4 * num_markers;
    size_t nbytes_marker_corrs = num_markers * (num_markers - 1) * 2;
    size_t nbytes_phen_data = num_individuals * num_phen * 4;
    size_t nbytes_marker_phen_corrs = num_markers * num_phen * 4;
    size_t req_mem_bytes =
        nbytes_marker_data + nbytes_marker_corrs + nbytes_phen_data + nbytes_marker_phen_corrs;
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    if (req_mem_bytes > device_mem_bytes)
    {
        // TODO: this is probably not exactly correct
        // figure out batch size
        float b = (float)num_individuals / 4.0;
        size_t max_batch_size =
            (size_t)(-b + std::sqrt(
                              b * b - (4.0 * (float)num_individuals * (float)num_phen -
                                       (float)device_mem_bytes)
                          ));
        size_t row_width = max_batch_size / 2;

        printf("Device mem < required mem; Running tiled routine. \n");

        cu_corr_pearson_npn_batched(
            marker_vals.data(),
            phen_vals.data(),
            num_markers,
            num_individuals,
            num_phen,
            marker_means.data(),
            marker_stds.data(),
            row_width,
            marker_corr.data(),
            marker_phen_corr.data(),
            phen_corr.data()
        );
        // TODO: the output corrs matrices are ordered in a really stupid way,
        // they need to be reordered.
    }
    else
    {
        printf("Calling correlation main\n");
        // compute correlations
        cu_corr_pearson_npn(
            marker_vals.data(),
            phen_vals.data(),
            num_markers,
            num_individuals,
            num_phen,
            marker_means.data(),
            marker_stds.data(),
            marker_corr.data(),
            marker_phen_corr.data(),
            phen_corr.data()
        );
    }

    printf("Writing results\n");
    // write results
    write_floats_to_binary(
        marker_corr.data(), marker_corr_mat_size, make_path(out_dir, chr_id, "_marker_corr.bin")
    );
    write_floats_to_binary(
        marker_phen_corr.data(),
        marker_phen_corr_mat_size,
        make_path(out_dir, chr_id, "_marker_phen_corr.bin")
    );
    write_floats_to_binary(
        phen_corr.data(), phen_corr_mat_size, make_path(out_dir, chr_id, "_phen_corr.bin")
    );

    std::cout << "Correlation matrices written to: " << out_dir << std::endl;
}

const std::string PRINTBF_USAGE = R"(
Interpret contents of binary file as floats and print to stoud.

usage: mps printbf <file>

arguments:
    file    file with floats in binary
)";

void printbf(int argc, char *argv[])
{
    // check for correct number of args
    if (argc != 3)
    {
        std::cout << PRINTBF_USAGE << std::endl;
        exit(1);
    }

    // check that path is valid
    std::string fpath = (std::string)argv[2];
    if (!path_exists(fpath))
    {
        std::cout << "file or directory found: " << fpath << std::endl;
        exit(1);
    }

    std::vector<float> read = read_floats_from_binary(fpath);

    // Using a for loop with index
    for (std::size_t i = 0; i < read.size(); ++i)
    {
        std::cout << read[i] << "\n";
    }
}

const std::string MPS_USAGE = R"(
usage: mps <command> [<args>]

commands:
    prep                    Prepare input (PLINK) .bed file for mps
    block                   Tile marker x marker correlation matrix
    cusk                    Run cuda-skeleton on block diagonal genomic covariance matrix
    cusk-single             Run cuda-skeleton on a single block of block diagonal genomic covariance matrix
    cuskss                  Run cuda-skeleton on a block of markers and traits with pre-computed correlations.
    cuskss-merged           Run cuda-skeleton on pre-selected markers and traits with pre-computed correlations.
    cuskss-trait-only       Run cuda-skeleton on a set of pre-computed trait-trait correlations.
    cusk-sim                Run cuda-skeleton on single simulated block
    cusk-phen               Run cuda-skeleton on phenotypes only

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
        cuda_skeleton_summary_stats(argc, argv);
    }
    else if (cmd == "cuskss-het")
    {
        cuda_skeleton_summary_stats_hetcor(argc, argv);
    }
    else if (cmd == "cuskss-merged")
    {
        cuda_skeleton_summary_stats_merged_blocks(argc, argv);
    }
    else if (cmd == "cuskss-trait-only")
    {
        cuskss_trait_only(argc, argv);
    }
    else if (cmd == "prep")
    {
        prep_bed(argc, argv);
    }
    else if (cmd == "cusk")
    {
        cusk(argc, argv);
    }
    else if (cmd == "cusk-single")
    {
        cusk_single(argc, argv);
    }
    else if (cmd == "block")
    {
        make_blocks(argc, argv);
    }
    else if (cmd == "cusk-phen")
    {
        phenotype_pc(argc, argv);
    }
    else if (cmd == "cusk-sim")
    {
        sim_pc(argc, argv);
    }
    else
    {
        std::cout << MPS_USAGE << std::endl;
    }

    return EXIT_SUCCESS;
}
