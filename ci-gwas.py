#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
from cusk_postprocessing.sepselect import sepselect_merged
from cusk_postprocessing.merge_blocks import (
    merge_block_outputs,
    reformat_cuskss_merged_output,
)

script_dir = os.path.dirname(os.path.realpath(__file__))
MPS_PATH = f"{script_dir}/cusk/build/apps/mps"
MVIVW_PATH = f"{script_dir}/mvivw/cig_mvivw.R"
RFCI_PATH = f"{script_dir}/srfci/CIGWAS_est_PAG.R"
DAVS_PATH = f"{script_dir}/sdavs/CIGWAS_est_ACE.R"


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


class TypeCheck:
    def __init__(self, type_fn, name: str, min_val=None, max_val=None):
        self._type_fn = type_fn
        self._name = name
        self._min_val = min_val
        self._max_val = max_val

    def __call__(self, val):
        val = self._type_fn(val)
        if self._min_val is not None and val < self._min_val:
            raise argparse.ArgumentTypeError(f"Minimum {self._name} is {self._min_val}")
        if self._max_val is not None and val > self._max_val:
            raise argparse.ArgumentTypeError(f"Maximum {self._name} is {self._max_val}")
        return val


def main():
    ci_gwas_parser = MyParser(
        prog="ci-gwas",
        description="Causal inference for multiple risk factors and diseases from genomics data",
    )
    subparsers = ci_gwas_parser.add_subparsers(
        required=True,
        title="subcommands",
        help="all individuals steps in the ci-gwas workflow",
    )

    # prep-bed
    prep_bed_parser = subparsers.add_parser(
        "prep-bed", help="Prepare PLINK bed file for cusk"
    )
    prep_bed_parser.add_argument(
        "bfiles", help="filestem of .bed, .bim, .fam fileset", type=str
    )
    prep_bed_parser.set_defaults(func=prep_bed)

    # block
    block_parser = subparsers.add_parser(
        "block",
        help="Tile whole-genome LD matrix into block diagonal matrix (requires GPU)",
    )
    block_parser.add_argument(
        "bfiles", help="filestem of .bed, .bim, .fam fileset", type=str
    )
    block_parser.add_argument(
        "max_block_size",
        metavar="max-block-size",
        help="maximum number of markers per block",
        default=11000,
        type=TypeCheck(int, "max-block-size", 2, None),
    )
    block_parser.add_argument(
        "device_mem_gb",
        metavar="device-mem-gb",
        help="maximum memory available on GPU in GB",
        default=10,
        type=TypeCheck(int, "device-mem-gb", 0, None),
    )
    block_parser.add_argument(
        "corr_width",
        metavar="corr-width",
        help="width of banded-correlation matrix",
        default=2000,
        type=TypeCheck(int, "corr-width", 2, None),
    )
    block_parser.set_defaults(func=block)

    # cusk (cusk-single)
    cusk_parser = subparsers.add_parser(
        "cusk",
        help="Infer skeleton with markers and traits as nodes, using marker data (requires GPU)",
    )
    cusk_parser.add_argument(
        "block_index",
        metavar="block-index",
        type=TypeCheck(int, "block-index", 0, None),
        help="0-based index of the block to run cusk on",
    )
    cusk_parser.add_argument(
        "blocks",
        help="file with genomic block definitions (output of ci-gwas block)",
        type=str,
    )
    cusk_parser.add_argument(
        "bfiles", help="filestem of .bed, .bim, .fam fileset", type=str
    )
    cusk_parser.add_argument(
        "phen", help="path to standardized phenotype tsv.", type=str
    )
    cusk_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    cusk_parser.add_argument(
        "max_level",
        metavar="max-level",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in cuPC (<= 14)",
        default=3,
    )
    cusk_parser.add_argument(
        "max_level_two",
        metavar="max-level-two",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in the second round of cuPC (<= 14)",
        default=14,
    )
    cusk_parser.add_argument(
        "max_depth",
        metavar="max-depth",
        type=TypeCheck(int, "max-depth", 1, None),
        help="max depth at which marker variables are kept as ancestors (>= 1)",
        default=1,
    )
    cusk_parser.add_argument(
        "outdir",
        type=str,
        help="directory for output",
        default="./",
    )
    cusk_parser.set_defaults(func=cusk)

    # cuskss
    cuskss_parser = subparsers.add_parser(
        "cuskss",
        help="Infer skeleton with markers and traits as nodes, using summary statistic data (requires GPU)",
    )
    cuskss_parser.add_argument(
        "mxm",
        type=str,
        help="Correlations between markers in block. Binary of floats, lower triangular, with diagonal, row major.",
    )
    cuskss_parser.add_argument(
        "mxp",
        type=str,
        help="Correlations between markers in all blocks and all traits. Textfile, whitespace separated, with columns: [chr, snp, ref, ...<trait names>], rectangular.",
    )
    cuskss_parser.add_argument(
        "pxp",
        type=str,
        help="Correlations between all traits. Textfile, whitespace separated, rectangular, only upper triangle is used. With trait names as column and row names. Order of traits has to be same as in the mxp file.",
    )
    cuskss_parser.add_argument(
        "block_index",
        metavar="block-index",
        type=TypeCheck(int, "block-index", 0, None),
        help="0-based index of the block to run cusk on",
    )
    cuskss_parser.add_argument(
        "blocks",
        help="file with genomic block definitions (output of ci-gwas block)",
        type=str,
    )
    cuskss_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    cuskss_parser.add_argument(
        "max_level",
        metavar="max-level",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in the first round of cuPC (<= 14)",
        default=3,
    )
    cuskss_parser.add_argument(
        "max_level_two",
        metavar="max-level-two",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in the second round of cuPC (<= 14)",
        default=14,
    )
    cuskss_parser.add_argument(
        "max_depth",
        metavar="max-depth",
        type=TypeCheck(int, "max-depth", 1, None),
        help="max depth at which marker variables are kept as ancestors (>= 1)",
        default=1,
    )
    cuskss_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    cuskss_parser.add_argument(
        "outdir",
        type=str,
        help="directory for output",
        default="./",
    )
    cuskss_parser.set_defaults(func=cuskss)

    # cuskss-merged
    cuskss_merged_parser = subparsers.add_parser(
        "cuskss-merged",
        help="Infer skeleton with markers and traits as nodes, using summary statistic data (requires GPU)",
    )
    cuskss_merged_parser.add_argument(
        "mxm",
        type=str,
        help="Correlations between selected markers. Binary of floats, lower triangular, with diagonal, row major.",
    )
    cuskss_merged_parser.add_argument(
        "mxp",
        type=str,
        help="Correlations between markers in all blocks and all traits. Textfile, whitespace separated, with columns: [chr, snp, ref, ...<trait names>], rectangular.",
    )
    cuskss_merged_parser.add_argument(
        "pxp",
        type=str,
        help="Correlations between all traits. Textfile, whitespace separated, rectangular, only upper triangle is used. With trait names as column and row names. Order of traits has to be same as in the mxp file.",
    )
    cuskss_merged_parser.add_argument(
        "marker_indices",
        metavar="marker-indices",
        type=str,
        help="Row indices if selected markers in mxp file. E.g. the .ixs file produced by `ci-gwas merge-block-outputs`. Binary of 32 bit ints.",
    )
    cuskss_merged_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    cuskss_merged_parser.add_argument(
        "max_level",
        metavar="max-level",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in the first round of cuPC (<= 14)",
        default=3,
    )
    cuskss_merged_parser.add_argument(
        "max_level_two",
        metavar="max-level-two",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in the second round of cuPC (<= 14)",
        default=14,
    )
    cuskss_merged_parser.add_argument(
        "max_depth",
        metavar="max-depth",
        type=TypeCheck(int, "max-depth", 1, None),
        help="max depth at which marker variables are kept as ancestors (>= 1)",
        default=1,
    )
    cuskss_merged_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    cuskss_merged_parser.add_argument(
        "outdir",
        type=str,
        help="directory for output",
        default="./",
    )
    cuskss_merged_parser.set_defaults(func=cuskss_merged)

    # merge block outputs
    merge_blocks_parser = subparsers.add_parser(
        "merge-block-outputs",
        help="Merge outputs from multiple blocks processed with cusk or cuskss",
    )
    merge_blocks_parser.add_argument(
        "cusk_output_dir",
        metavar="cusk-output-dir",
        type=str,
        help="output directory of cusk or cuskss",
    )
    merge_blocks_parser.add_argument(
        "blockfile",
        help="file with genomic block definitions (output of ci-gwas block)",
        type=str,
    )
    merge_blocks_parser.set_defaults(func=merge_blocks)

    # mvivw
    mvivw_parser = subparsers.add_parser(
        "mvivw",
        help="Run multivariable inverse-variance weighted mendelian randomization between all adjacent traits, using cusk-identified markers as intrumental variables",
    )
    mvivw_parser.add_argument(
        "cusk_output_dir",
        metavar="cusk-output-dir",
        type=str,
        help="output directory of cusk or cuskss",
    )
    mvivw_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    mvivw_parser.set_defaults(func=run_mvivw)

    # sepselect
    sepselect_parser = subparsers.add_parser(
        "sepselect",
        help="Compute maximal and partial-correlation-minimizing separation sets on merged cusk skeletons",
    )
    sepselect_parser.add_argument(
        "cusk_result_stem",
        metavar="cusk-result-stem",
        help="outdir + stem of merged cusk results",
        type=str,
    )
    sepselect_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    sepselect_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    sepselect_parser.set_defaults(func=run_sepselect)

    # sRFCI
    srfci_parser = subparsers.add_parser("srfci", help="Run sRFCI to infer a PAG")
    srfci_parser.add_argument(
        "sepselect_result_stem",
        metavar="sepselect-result-stem",
        help="outdir + stem of sepselect",
        type=str,
    )
    srfci_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    srfci_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    srfci_parser.set_defaults(func=srfci)

    # sDAVS
    sdavs_parser = subparsers.add_parser(
        "sdavs", help="Run sDAVS on a pair of traits to infer the ACE"
    )
    sdavs_parser.add_argument(
        "exposure",
        help="1-based index of exposure trait",
        type=TypeCheck(int, "exposure", 1, None),
    )
    sdavs_parser.add_argument(
        "outcome",
        help="1-based index of outcome trait",
        type=TypeCheck(int, "outcome", 1, None),
    )
    sdavs_parser.add_argument(
        "pag_path",
        metavar="pag-path",
        help="path to pag estimated with srfci",
        type=str,
    )
    sdavs_parser.add_argument(
        "output_file",
        metavar="output-file",
        help="output filename",
        type=str,
    )
    sdavs_parser.add_argument(
        "sepselect_result_stem",
        metavar="sepselect-result-stem",
        help="outdir + stem of sepselect",
        type=str,
    )
    sdavs_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    sdavs_parser.add_argument(
        "num_samples",
        metavar="num-samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    sdavs_parser.set_defaults(func=sdavs)

    args = ci_gwas_parser.parse_args()
    args.func(args)


def prep_bed(args):
    subprocess.run([MPS_PATH, "prep", args.bfiles], check=True)


def block(args):
    subprocess.run(
        [
            MPS_PATH,
            "block",
            args.bfiles,
            str(args.max_block_size),
            str(args.device_mem_gb),
            str(args.corr_width),
        ],
        check=True,
    )


def cusk(args):
    # check that everything is standardized, and that all files exist
    subprocess.run(
        [
            MPS_PATH,
            "cusk-single",
            args.phen,
            args.bfiles,
            args.blocks,
            str(args.alpha),
            str(args.max_level),
            str(args.max_level_two),
            str(args.max_depth),
            args.outdir,
            str(args.block_index),
        ],
        check=True,
    )


def cuskss(args):
    subprocess.run(
        [
            MPS_PATH,
            "cuskss",
            args.mxm,
            args.mxp,
            args.pxp,
            str(args.block_index),
            args.blocks,
            str(args.alpha),
            str(args.max_level),
            str(args.max_level_two),
            str(args.max_depth),
            str(args.num_samples),
            args.outdir,
        ],
        check=True,
    )


def cuskss_merged(args):
    subprocess.run(
        [
            MPS_PATH,
            "cuskss-merged",
            args.mxm,
            args.mxp,
            args.pxp,
            args.marker_indices,
            str(args.alpha),
            str(args.max_level),
            str(args.max_level_two),
            str(args.max_depth),
            str(args.num_samples),
            args.outdir,
        ],
        check=True,
    )
    # reformat the output to conform with the merge_blocks format
    reformat_cuskss_merged_output(cusk_dir=args.outdir).write_mm(
        basepath=f"{args.outdir}/cuskss_merged"
    )


def merge_blocks(args):
    out_dir = args.cusk_output_dir
    if not out_dir.endswith("/"):
        out_dir += "/"
    merged_blocks = merge_block_outputs(args.blockfile, out_dir)
    merged_blocks.write_mm(f"{args.cusk_output_dir}/merged_blocks")


def run_sepselect(args):
    merged_cusk = sepselect_merged(args.cusk_result_stem, args.alpha, args.num_samples)
    merged_cusk.to_file(f"{os.path.dirname(args.cusk_result_stem)}/max_sep_min_pc")


def run_mvivw(args):
    subprocess.run(
        [
            MVIVW_PATH,
            args.cusk_output_dir,
            str(args.num_samples),
            "TRUE", # rm exposures that have been identified as non-adjacent in cusk
            "FALSE", # use ld matrix
            f"{args.cusk_output_dir}/mvivw_results.tsv", # output dir
        ],
        check=True,
    )


def srfci(args):
    subprocess.run(
        [
            "Rscript",
            RFCI_PATH,
            args.sepselect_result_stem,
            str(args.alpha),
            str(args.num_samples),
            "cusk2",
        ],
        check=True,
    )


def sdavs(args):
    subprocess.run(
        [
            "Rscript",
            DAVS_PATH,
            str(args.exposure),
            str(args.outcome),
            str(args.num_samples),
            str(args.alpha),
            args.pag_path,
            args.sepselect_result_stem,
            args.output_file,
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
