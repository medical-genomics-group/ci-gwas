auto main(int argc, char *argv[]) -> int
{
    if (argc == 1)
    {
        return EXIT_SUCCESS;
    }
}

// #include <mps/cli.h>
// #include <mps/io.h>

// const std::string MPS_USAGE = R"(
// usage: mps <command> [<args>]

// commands:
//     prep                    Prepare input (PLINK) .bed file for mps
//     block                   Tile marker x marker correlation matrix
//     cusk                    Run cuda-skeleton on a single block of block diagonal genomic covariance matrix
//     cuskss                  Run cuda-skeleton on a block of markers and traits with pre-computed correlations.

// contact:
//     nick.machnik@gmail.com
// )";

// auto main(int argc, char *argv[]) -> int
// {
//     if (argc == 1)
//     {
//         std::cout << MPS_USAGE << std::endl;
//         return EXIT_SUCCESS;
//     }

//     std::string cmd = (std::string)argv[1];

//     if ((cmd == "--help") || (cmd == "-h"))
//     {
//         std::cout << MPS_USAGE << std::endl;
//     }
//     else if (cmd == "cuskss")
//     {
//         std::string mxm_path = argv[2];
//         std::string mxp_path = argv[3];
//         std::string mxp_se_path = argv[4];
//         std::string pxp_path = argv[5];
//         std::string pxp_se_path = argv[6];
//         std::string time_index_path = argv[7];
//         int block_index = std::stoi(argv[8]);
//         std::string blockfile_path = argv[9];
//         std::string marker_index_path = argv[10];
//         float alpha = std::stof(argv[11]);
//         int max_level_one = std::stoi(argv[12]);
//         int max_level_two = std::stoi(argv[13]);
//         int max_depth = std::stoi(argv[14]);
//         float num_samples = (float)std::stoi(argv[15]);
//         std::string outdir_path = argv[16];

//         bool merged = marker_index_path != "NULL";
//         bool hetcor = mxp_se_path != "NULL";
//         bool trait_only = mxm_path != "NULL";
//         bool two_stage = max_level_two > 0;
//         bool time_indexed = time_index_path != "NULL";

//         // required args
//         check_path(pxp_path);
//         check_path(blockfile_path);
//         check_path(outdir_path);

//         if (hetcor || merged || !trait_only) {
//             check_path(mxm_path);
//             check_path(mxp_path);
//         }

//         if (hetcor) {
//             check_path(mxp_se_path);
//             check_path(pxp_se_path);
//         }

//         if (time_indexed) {
//             check_path(time_index_path);
//         }

//         CuskssArgs args = {
//             merged,             // merged
//             hetcor,             // hetcor
//             trait_only,         // trait_only
//             two_stage,          // two_stage
//             time_indexed,       // time_indexed
//             alpha,              // alpha
//             num_samples,        // pearson_sample_size
//             max_level_one,      // max_level
//             max_level_two,      // max_level_two
//             max_depth,          // depth
//             block_index,        // block_ix
//             blockfile_path,         // block_path
//             marker_index_path,  // marker_ixs_path
//             mxm_path,           // mxm_path
//             mxp_path,           // mxp_path
//             mxp_se_path,        // mxp_se_path
//             pxp_path,           // pxp_path
//             pxp_se_path,        // pxp_se_path
//             time_index_path,    // time_index_path
//             outdir_path         // outdir
//         };
//         cuskss(args);
//     }
//     else if (cmd == "prep")
//     {
//         prep_bed(argc, argv);
//     }
//     else if (cmd == "cusk")
//     {
//         cusk(argc, argv);
//     }
//     else if (cmd == "block")
//     {
//         make_blocks(argc, argv);
//     }
//     else
//     {
//         std::cout << MPS_USAGE << std::endl;
//     }

//     return EXIT_SUCCESS;
// }
