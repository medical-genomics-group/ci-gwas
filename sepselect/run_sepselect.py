#!/usr/bin/env python

import sys
import sepselect as sepselect


def main():
    results_dir = sys.argv[1]
    num_samples = int(sys.argv[2])
    alpha_e = int(sys.argv[3])
    cr = sepselect.sepselect_merged(
        results_dir + "/all_merged", 10.0 ** (-alpha_e), num_samples
    )
    cr.to_file(results_dir + "/max_sep_min_pc")


if __name__ == "__main__":
    main()
