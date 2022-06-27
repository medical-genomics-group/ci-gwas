#include <iostream>
#include <string>
#include <cassert>
#include <sys/stat.h>

#include <mps/prep_markers.h>

const std::string MPS_USAGE = R"(
usage: mps <command> [<args>]

commands:
    prep    Prepare input (PLINK) .bed file for mps
    corr    Compute the marker/phenotype correlation matrix
    cups    use cuPC to compute the parent set for each phenotype

contact:
    nick.machnik@gmail.com
)";

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
        struct stat buffer;
        if (stat (((std::string)argv[i]).c_str(), &buffer) != 0) {
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

