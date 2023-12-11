# CI-GWAS: Causal inference for multiple risk factors and diseases from genomics data

This is the official implementation of CI-GWAS as described in our paper: [Causal inference for multiple risk factors and diseases from genomics data](https://www.biorxiv.org/content/10.1101/2023.12.06.570392v1)

## Prerequisites

- `cmake >= 3.18`
- `python >= 3.9.6` , with `numpy >= 1.22.1` and `scipy >= 1.11.3`.
- `R >= 4.1.0`, with dependencies:

```R
install.packages("BiocManager")
BiocManager::install(pkgs=c("graph","Rgraphviz", "RBGL"))
install.packages(c( "abind", "igraph", "ggm", "corpcor", "robustbase", "vcd", "Rcpp", "bdsmatrix", "sfsmisc", "fastICA", "clue", "MASS", "Matrix", "mvtnorm"," huge", "ggplot2", "dagitty", 'pcalg', 'Matrix'))
```

## Installation

CI-GWAS is a conglomerate of scripts and compiled programs bundled in a python command-line-interface.

First simply clone the repo:

```bash
git clone --recurse-submodules https://github.com/medical-genomics-group/ci-gwas.git
```

The cli should already be accessible via

```bash
./ci-gwas.py
```

The `cusk` part of the project has to compiled:

```bash
cd cusk
cmake -S . -B build
cmake --build build
```

You can then run the tests to check that everything works:

```bash
cd build && ctest
```

## Running

Get help:

```bash
./ci-gwas.py -h
```

First of all, make sure that any marker data you want to plug in is LD pruned, or at least does not have markers with a correlation of 1.

A standard analysis, if you have data at the individual-level available, consists of subsquent calls to

1. `ci-gwas.py prep-bed` to compute means and variances of all markers
2. `ci-gwas.py block` to block the LD matrix
3. `ci-gwas.py cusk` (once for each block) to compute skeletons
   (**make sure that the trait values are standardized**)
4. `ci-gwas.py merge-block-outputs` to merge all skeletons
5. `ci-gwas.py sepselect` to find separation sets
6. `ci-gwas.py srfci` to infer a PAG
7. `ci-gwas.py srfci` to infer ACEs

Alternatively, if you have correlations from summarized data, you can start at step 3) with `cuskss` instead of `cusk`. In that case it is important that

- the traits have the same order in the `mxp` and `pxp` files
- the markers have the same order in the `mxm` and `mxp` files

## Common Errors

`no kernel image is available for execution on the device`

Can be caused by the chosen device having a lower [GPU Compute Capability](https://developer.nvidia.com/cuda-gpus#collapseOne) than the one `cusk` was compiled for. The Compute Capability targeted by the build is specified in the top level `CMakeLists.txt`. If there are multiple devices on the machine and only a subset of them have the appropriate Compute Capability, you can choose one by setting `CUDA_VISIBLE_DEVICES=X` where X is the index of the device. The device list can be inspected with `nvidia-smi -L`.
