#include <iostream>
#include <string>
#include <cassert>
#include <sys/stat.h>

#include <mps/corr_compressed.h>
#include <mps/prep_markers.h>
#include <mps/io.h>

const std::string PREP_USAGE = R"(
Prepare input (PLINK) .bed file for mps.

The file is split by chromosome, NaNs are imputed to column medians
and column means and standard deviations are computed.

usage: mps prep <.bed> <.bim> <.fam> <outdir> <mem_gb>

arguments:
    .bed    (PLINK) .bed file with marker data
    .bim    (PLINK) .bim file with allele information
    .fam    (PLINK) .fam file with sample information
    outdir  output directory for processed .bed
    mem_gb  maximal amount of memory available in Gb
)";

auto path_exists(std::string path) -> bool
{
    struct stat buffer;
    return (stat (path.c_str(), &buffer) == 0);
}

void prep_bed(int argc, char *argv[])
{
    // check if enough args present
    if ((argc != 7) || (argv[2] == "--help") || (argv[2] == "-h")) 
    {
        std::cout << PREP_USAGE << std::endl;
        exit(1);
    }

    // check if files and dirs exist
    for (size_t i = 2; i < 6; ++i) 
    {
        if (! path_exists((std::string)argv[i])) {
            std::cout << "file or directory found: " << argv[i] << std::endl;
            exit(1);
        }
    }

    int mem_gb = stoi((std::string)argv[6]);
    std::string bed_path = (std::string)argv[2];
    std::string bim_path = (std::string)argv[3];
    std::string fam_path = (std::string)argv[4];
    std::string out_dir = (std::string)argv[5];

    // run prep
    prep_bed(bed_path, bim_path, fam_path, out_dir, mem_gb);

    std::cout << "Preprocessing output files written to: " << out_dir << std::endl;
}

const std::string CORR_USAGE = R"(
Compute correlations between markers and phenotypes.

usage: mps corr <prepdir> <chr> <device_mem_gb> <.phen>...

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dims files for each chromosome
    chr ID of the chromsome to be processed
    device_mem_gb   Amount of memory available on the GPU
    .phen   Path to .phen file with phenotype values, sorted in the same way as genotype info in the original bed file
)";

void corr(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < 6) {
        std::cout << CORR_USAGE << std::endl;
        exit(1);
    }

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dims", ".means", ".stds"};
    for (size_t i = 0; i < 4; ++i) {
        std::string fpath = make_path(out_dir, chr_id, req_suffixes[i]);
        if (!path_exists(fpath)) {
            std::cout << "file or directory found: " << fpath << std::endl;
            exit(1);
        }
    }

    std::vector<std::string> phen_paths = {};
    for (size_t i = 5; i < argc; ++i) {
        std::string phen_path = (std::string)argv[i];
        if (!path_exists(phen_path)) {
            std::cout << "file or directory found: " << phen_path << std::endl;
            exit(1);
        }
        phen_paths.push_back(phen_path);
    }
    
    // load data + check that file contents are valid
    // .dims
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dims");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2) {
        std::cout << "Invalid .dims file: found " << ndims << "dimensions instead of two." <<  std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    // .phen
    size_t num_phen = phen_paths.size();
    std::vector<float> phen = {};
    for (size_t i = 0; i < num_phen; ++i) {
        read_floats_from_lines(phen_paths[i], phen);
    }

    // .bed
    std::string bed_path = make_path(out_dir, chr_id, ".bed");
    size_t nbytes_per_block = (num_individuals + 3) / 4;
    size_t nbytes_marker_vals = num_markers * nbytes_per_block;
    std::vector<unsigned char> marker_vals(nbytes_marker_vals);
    read_n_bytes_from_binary(bed_path, nbytes_marker_vals, marker_vals);

    //size_t num_individuals = dims[0];
    //size_t num_markers = dims[1];

    //const size_t num_markers = BMT_NUM_MARKERS;
    //const size_t num_individuals = BMT_NUM_INDIVIDUALS;
    //const size_t num_phen = BMT_NUM_PHEN;
    //const size_t marker_cm_size = corr_matrix_size(num_markers);
    //const size_t marker_phen_cm_size = num_markers * num_phen;
    //const size_t phen_cm_size = corr_matrix_size(num_phen);
 
    //float marker_corr[marker_cm_size];
    //memset(marker_corr, 0.0, sizeof(marker_corr));
    //float marker_phen_corr[marker_phen_cm_size];
    //memset(marker_phen_corr, 0.0, sizeof(marker_phen_corr));
    //float phen_corr[phen_cm_size];
    //memset(phen_corr, 0.0, sizeof(phen_corr));

    //cu_corr_npn(bmt_marker_vals, bmt_phen_vals, num_markers, num_individuals, num_phen,
    //             bmt_marker_mean, bmt_marker_std, marker_corr, marker_phen_corr, phen_corr);
}

const std::string MPS_USAGE = R"(
usage: mps <command> [<args>]

commands:
    prep    Prepare input (PLINK) .bed file for mps
    corr    Compute the marker/phenotype correlation matrix
    cups    use cuPC to compute the parent set for each phenotype

contact:
    nick.machnik@gmail.com
)";

auto main(int argc, char *argv[]) -> int 
{
    if (argc == 1) {
        std::cout << MPS_USAGE << std::endl;
        return EXIT_SUCCESS;
    }
    
    std::string cmd = (std::string)argv[1];
    
    if ((cmd == "--help") || (cmd == "-h")) {
        std::cout << MPS_USAGE << std::endl;
    } else if (cmd == "prep") {
        prep_bed(argc, argv);
    } else if (cmd == "corr") {
        std::cout << "'corr' cli is not implemented yet." << std::endl;
    } else if (cmd == "cups") {
        std::cout << "'cups' cli is not implemented yet." << std::endl;
    } else {
        std::cout << MPS_USAGE << std::endl;
    }

    return EXIT_SUCCESS;
}

