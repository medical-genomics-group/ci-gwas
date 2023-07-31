#!/usr/bin/env python

import argparse
import sys
import subprocess


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
    ci_gwas_parser = argparse.ArgumentParser(
        prog="ci-gwas",
        description="Causal inference for multiple risk factors and diseases from summary statistic data",
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
        "max-block-size",
        help="maximum number of markers per block",
        default=11000,
        type=TypeCheck(int, "max-block-size", 2, None),
    )
    block_parser.add_argument(
        "device-mem-gb",
        help="maximum memory available on GPU in GB",
        default=10,
        type=TypeCheck(int, "device-mem-gb", 0, None),
    )
    block_parser.add_argument(
        "corr-width",
        help="width of banded-correlation matrix",
        default=2000,
        type=TypeCheck(int, "corr-width", 2, None),
    )
    block_parser.set_defaults(func=block)

    # cusk
    cusk_parser = subparsers.add_parser(
        "cusk",
        help="Infer skeleton with markers and traits as nodes, using marker data (requires GPU)",
    )
    cusk_parser.add_argument(
        "block-index",
        type=TypeCheck(int, "block-index", 0, None),
        type=str,
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
        "phen", help="path to standardized phenotype tsv", type=str
    )
    cusk_parser.add_argument(
        "alpha",
        type=TypeCheck(float, "alpha", 0.0, 1.0),
        help="significance level for conditional independence tests",
        default=10**-4,
    )
    cusk_parser.add_argument(
        "max-level",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in cuPC (<= 14)",
        default=6,
    )
    cusk_parser.add_argument(
        "max-depth",
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
        help="Correlations between markers in block. Binary of floats, upper triangular, without diagonal.",
    )
    cuskss_parser.add_argument(
        "mxp",
        type=str,
        help="Correlations between markers in all blocks and all traits. Textfile, whitespace separated, with columns: [chr, snp, ref, ...<trait names>], rectangular.",
    )
    cuskss_parser.add_argument(
        "pxp",
        type=str,
        help="Correlations between all traits. Textfile, whitespace separated, rectangular, only upper triangle is used. With trait names as column and row names.",
    )
    cuskss_parser.add_argument(
        "block-index",
        type=TypeCheck(int, "block-index", 0, None),
        type=str,
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
        "max-level",
        type=TypeCheck(int, "max-level", 1, 14),
        help="maximal size of separation sets in cuPC (<= 14)",
        default=6,
    )
    cuskss_parser.add_argument(
        "max-depth",
        type=TypeCheck(int, "max-depth", 1, None),
        help="max depth at which marker variables are kept as ancestors (>= 1)",
        default=1,
    )
    cuskss_parser.add_argument(
        "num_samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
    )
    cuskss_parser.add_argument(
        "num_samples",
        type=TypeCheck(int, "num-samples", 1, None),
        help="number of samples used for computing correlations",
        default=1,
    )
    cuskss_parser.add_argument(
        "outdir",
        type=str,
        help="directory for output",
        default="./",
    )
    cuskss_parser.set_defaults(func=cuskss)

    # sRFCI
    srfci_parser = subparsers.add_parser("srfci", help="Run sRFCI to infer a PAG")
    srfci_parser.set_defaults(func=srfci)

    # sDAVS
    sdavs_parser = subparsers.add_parser("sdavs", help="Run sDAVS to infer ACEs")
    sdavs_parser.set_defaults(func=sdavs)

    args = ci_gwas_parser.parse_args(sys.argv[1:])
    args.func(args)


def prep_bed(args):
    res = subprocess.run(["mps", args.bfiles], check=True)


def block(args):
    res = subprocess.run(
        ["mps", args.bfiles, args.max_block_size, args.device_mem_gb, args.corr_width],
        check=True,
    )


def cusk(args):
    pass


def cuskss(args):
    pass


def srfci(args):
    pass


def sdavs(args):
    pass


if __name__ == "__main__":
    main()
