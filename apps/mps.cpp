#include <math.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/io.h>
#include <mps/prep_markers.h>
#include <sys/stat.h>

#include <array>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

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

auto path_exists(const std::string path) -> bool
{
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// TODO: this should also parse .phen files and make sure that no info is missing and order checks
// out
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
        if (!path_exists((std::string)argv[i]))
        {
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

const std::string MCORRK_USAGE = R"(
Compute Pearson correlation between markers as sin(pi / 2 * tau-b) where tau is the Kendall's tau-b

usage: mps mcorrk <prepdir> <chr> <device_mem_gb>

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dims files for each chromosome
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

    double device_mem_gb = atof(argv[4]);
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dims"};
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
    // .dims
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dims");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dims file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

    printf("Loading .bed \n");

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

    if (req_mem_bytes < device_mem_bytes)
    {
        // figure out batch size
        double b = num_individuals / 4 - 2;
        size_t max_batch_size = (size_t)((-b + std::sqrt(b * b + 8 * device_mem_bytes)) / 4);
        size_t batch_nrows = max_batch_size / 2;

        printf("Device mem < required mem; Running tiled routine. \n");

        cu_marker_corr_pearson_npn_batched(marker_vals.data(), num_markers, num_individuals,
                                           batch_nrows, marker_corr.data());
    }
    else
    {
        printf("Calling correlation main\n");
        // compute correlations
        cu_marker_corr_pearson_npn(marker_vals.data(), num_markers, num_individuals,
                                   marker_corr.data());
    }

    printf("Writing results\n");
    // write results
    write_floats_to_binary(marker_corr.data(), marker_corr_mat_size,
                           make_path(out_dir, chr_id, "_marker_corrk.bin"));

    std::cout << "Correlation matrix written to: " << out_dir << std::endl;
}

const std::string MCORRP_USAGE = R"(
Compute Pearson correlation between markers

usage: mps mcorrp <prepdir> <chr> <device_mem_gb>

arguments:
    prepdir Directory with `mps prep` output, i.e. .bed, .stds, .means and .dims files for each chromosome
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

    double device_mem_gb = atof(argv[4]);
    size_t device_mem_bytes = device_mem_gb * std::pow(10, 9);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dims", ".means", ".stds"};
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
    // .dims
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dims");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dims file: found " << ndims << "dimensions instead of two."
                  << std::endl;
        exit(1);
    }
    size_t num_individuals = dims[0];
    size_t num_markers = dims[1];

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

    // allocate correlation result arrays
    size_t marker_corr_mat_size = num_markers * (num_markers - 1) / 2;
    std::vector<float> marker_corr(marker_corr_mat_size, 0.0);

    // mem required for non-batched processing
    size_t nbytes_marker_data = num_individuals / 4 * num_markers;
    size_t nbytes_marker_corrs = num_markers * (num_markers - 1) * 2;
    size_t req_mem_bytes = nbytes_marker_data + nbytes_marker_corrs;

    if (req_mem_bytes < device_mem_bytes)
    {
        // figure out batch size
        double b = num_individuals / 4 - 2;
        size_t max_batch_size = (size_t)((-b + std::sqrt(b * b + 8 * device_mem_bytes)) / 4);
        size_t batch_nrows = max_batch_size / 2;

        printf("Device mem < required mem; Running tiled routine. \n");

        cu_marker_corr_pearson_batched(marker_vals.data(), num_markers, num_individuals,
                                       marker_means.data(), marker_stds.data(), batch_nrows,
                                       marker_corr.data());
    }
    else
    {
        printf("Calling correlation main\n");
        // compute correlations
        cu_marker_corr_pearson(marker_vals.data(), num_markers, num_individuals,
                               marker_means.data(), marker_stds.data(), marker_corr.data());
    }

    printf("Writing results\n");
    // write results
    write_floats_to_binary(marker_corr.data(), marker_corr_mat_size,
                           make_path(out_dir, chr_id, "_marker_corrp.bin"));

    std::cout << "Correlation matrix written to: " << out_dir << std::endl;
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

const int CORR_NARGS = 6;

void corr(int argc, char *argv[])
{
    // check for correct number of args
    if (argc < CORR_NARGS)
    {
        std::cout << CORR_USAGE << std::endl;
        exit(1);
    }

    double device_mem_gb = atof(argv[4]);

    // check that paths are valid
    // TODO: figure out how to glob files and check that at least one .phen is present
    std::string out_dir = (std::string)argv[2];
    std::string chr_id = (std::string)argv[3];
    std::string req_suffixes[4] = {".bed", ".dims", ".means", ".stds"};
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
    // .dims
    // TODO: there is a lot more that can go wrong here, e.g. number of cols
    std::string dims_path = make_path(out_dir, chr_id, ".dims");
    std::vector<int> dims = read_ints_from_lines(dims_path);
    size_t ndims = dims.size();
    if (ndims != 2)
    {
        std::cout << "Invalid .dims file: found " << ndims << "dimensions instead of two."
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

    if (req_mem_bytes < device_mem_bytes)
    {
        // TODO: this is probably not exactly correct
        // figure out batch size
        double b = (double)num_individuals / 4.0;
        size_t max_batch_size =
            (size_t)(-b + std::sqrt(b * b - (4.0 * (double)num_individuals * (double)num_phen -
                                             (double)device_mem_bytes)));
        size_t row_width = max_batch_size / 2;

        printf("Device mem < required mem; Running tiled routine. \n");

        cu_corr_pearson_npn_batched(marker_vals.data(), phen_vals.data(), num_markers,
                                    num_individuals, num_phen, marker_means.data(),
                                    marker_stds.data(), row_width, marker_corr.data(),
                                    marker_phen_corr.data(), phen_corr.data());
        // TODO: the output corrs matrices are ordered in a really stupid way,
        // they need to be reordered.
    }
    else
    {
        printf("Calling correlation main\n");
        // compute correlations
        cu_corr_pearson_npn(marker_vals.data(), phen_vals.data(), num_markers, num_individuals,
                            num_phen, marker_means.data(), marker_stds.data(), marker_corr.data(),
                            marker_phen_corr.data(), phen_corr.data());
    }

    printf("Writing results\n");
    // write results
    write_floats_to_binary(marker_corr.data(), marker_corr_mat_size,
                           make_path(out_dir, chr_id, "_marker_corr.bin"));
    write_floats_to_binary(marker_phen_corr.data(), marker_phen_corr_mat_size,
                           make_path(out_dir, chr_id, "_marker_phen_corr.bin"));
    write_floats_to_binary(phen_corr.data(), phen_corr_mat_size,
                           make_path(out_dir, chr_id, "_phen_corr.bin"));

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
    prep    Prepare input (PLINK) .bed file for mps
    corr    Compute the marker/phenotype correlation matrix
    mcorrk  Compute pearson correlations between markers as sin(pi / 2 tau_b)
    mcorrp  Compute pearson correlations between markers
    cups    use cuPC to compute the parent set for each phenotype

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
    else if (cmd == "prep")
    {
        prep_bed(argc, argv);
    }
    else if (cmd == "corr")
    {
        corr(argc, argv);
    }
    else if (cmd == "cups")
    {
        std::cout << "'cups' cli is not implemented yet." << std::endl;
    }
    else if (cmd == "printbf")
    {
        printbf(argc, argv);
    }
    else if (cmd == "mcorrp")
    {
        mcorrp(argc, argv);
    }
    else if (cmd == "mcorrk")
    {
        mcorrk(argc, argv);
    }
    else
    {
        std::cout << MPS_USAGE << std::endl;
    }

    return EXIT_SUCCESS;
}
