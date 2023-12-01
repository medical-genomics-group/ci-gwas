# CI-GWAS: Causal inference for multiple risk factors and diseases from genomics data
This is the official implementation of CI-GWAS as described in our paper: [Causal inference for multiple risk factors and diseases from genomics data]()

## Installation

CI-GWAS is a conglomerate of scripts and compiled programs to be executed on the command line

### Prerequisites

For R 4.1.0:

```R
install.packages("BiocManager")
BiocManager::install(pkgs=c("graph","Rgraphviz", "RBGL"))
install.packages(c( "abind", "igraph", "ggm", "corpcor", "robustbase", "vcd", "Rcpp", "bdsmatrix", "sfsmisc", "fastICA", "clue", "MASS", "Matrix", "mvtnorm"," huge", "ggplot2", "dagitty", 'pcalg', 'Matrix'))
```

## Running
