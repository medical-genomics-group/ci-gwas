#!/usr/bin/env python

import argparse
import sys
import subprocess


def main():
    ci_gwas_parser = argparse.ArgumentParser(
        prog="ci-gwas",
        description=
        "Causal inference for multiple risk factors and diseases from summary statistic data"
    )
    subparsers = ci_gwas_parser.add_subparsers(
        required=True,
        title="subcommands",
        help="all individuals steps in the ci-gwas workflow")

    # prep-bed
    prep_bed_parser = subparsers.add_parser(
        "prep-bed", help="Prepare PLINK bed file for cusk")
    prep_bed_parser.add_argument("bfiles",
                                 help="filestem of .bed, .bim, .fam fileset")
    prep_bed_parser.set_defaults(func=prep_bed)

    # block
    block_parser = subparsers.add_parser(
        "block",
        help=
        "Tile whole-genome LD matrix into block diagonal matrix (requires GPU)"
    )
    block_parser.add_argument("bfiles",
                              help="filestem of .bed, .bim, .fam fileset")
    block_parser.add_argument("max-block-size",
                              help="maximum number of markers per block")
    block_parser.add_argument("device-mem-gb",
                              help="maximum memory available on GPU in GB")
    block_parser.add_argument(
        "corr-width", help="max distance at which to compute correlations")
    block_parser.set_defaults(func=block)

    # cusk
    cusk_parser = subparsers.add_parser(
        "cusk",
        help=
        "Infer the skeleton with markers and traits as nodes, using marker data (requires GPU)"
    )
    cusk_parser.set_defaults(func=cusk)

    # cuskss
    cuskss_parser = subparsers.add_parser(
        "cuskss",
        help=
        "Infer the skeleton with markers and traits as nodes, using summary statistic data (requires GPU)"
    )
    cuskss_parser.set_defaults(func=cuskss)

    # sRFCI
    srfci_parser = subparsers.add_parser("srfci",
                                         help="Run sRFCI to infer a PAG")
    srfci_parser.set_defaults(func=srfci)

    # sDAVS
    sdavs_parser = subparsers.add_parser("sdavs",
                                         help="Run sDAVS to infer ACEs")
    sdavs_parser.set_defaults(func=sdavs)

    args = ci_gwas_parser.parse_args(sys.argv[1:])
    args.func(args)


def prep_bed(args):
    res = subprocess.run(["mps", args.bfiles], check=True)


def block(args):
    res = subprocess.run([
        "mps", args.bfiles, args.max_block_size, args.device_mem_gb,
        args.corr_width
    ],
                         check=True)


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
