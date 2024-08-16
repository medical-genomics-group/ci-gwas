library("pcalg")
library("fastmatch")
library("Matrix")

source("/nfs/scistore17/robingrp/nmachnik/dev/ci-gwas/srfci/RFCI_functions_final.R")

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)
input_filestem <- myargs[1]
alpha <- as.numeric(myargs[2])
num_individuals <- as.numeric(myargs[3])

print(paste0("input_filestem: ", input_filestem))
print(paste0("num_individuals: ", num_individuals))
print(paste0("alpha: ", alpha))

unsh_triple_pheno_only <- FALSE
force_marker_to_trait_before_R1 <- FALSE
force_marker_to_trait_in_the_end <- TRUE
external_ambiguous_triples <- TRUE
external_unshielded_triples <- TRUE

load_sparse_sepsets <- function(filepath, nvar) {
    res <- rep(list(rep(list(NULL), times = nvar)), nvar)
    con <- file(filepath, "r")
    while (TRUE) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) {
            break
        }
        vals <- as.integer(strsplit(line, " +")[[1]])
        gvals <<- vals
        res[[vals[1]]][[vals[2]]] <- tail(vals, -2)
    }
    close(con)
    return(res)
}

# number of variables, number of phenotypes, max_level
mdim_path <- paste0(input_filestem, ".mdim")
# separation set
sep_path <- paste0(input_filestem, ".ssm")
atr_path <- paste0(input_filestem, ".atr")
ut_path <- paste0(input_filestem, ".ut")

x <- read.csv(mdim_path, sep = "\t", header = FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
max_level <- as.numeric(x[3])
num_atr <- as.numeric(x[4])
num_ut <- as.numeric(x[5])
sepset <- load_sparse_sepsets(sep_path, num_var)

adjmat <- as.matrix(readMM(paste0(input_filestem, "_sam.mtx")))
cormat <- as.matrix(readMM(paste0(input_filestem, "_scm.mtx")))
A <- as.matrix(readMM(paste0(input_filestem, "_spm.mtx")))

atr_file <- file(atr_path, "rb")
atr <- readBin(atr_file, integer(), size = 4, n = num_atr * 3)
atrmat <- matrix(atr, ncol = 3, byrow = TRUE)

# mark unshielded triples a-b-c as unfaithful / ambiguous,
# if b is in maximal sepset, but not in minimal sepset.
p <- num_var
unfTrip = c()
for (i in 1:nrow(atrmat)) {
    # increment because of one based indexing
    x <- atrmat[i, 1] + 1
    y <- atrmat[i, 2] + 1
    z <- atrmat[i, 3] + 1
    c(unfTrip, triple2numb(p, x, y, z))
}

print("Applying R1-R10")
res <- udag2apag_ci_gwas(A[1:num_phen,1:num_phen], sepset, rules=rep(TRUE, 10), unfVect=unfTripl)
tmp_graph <- A
tmp_graph[1:num_phen, 1:num_phen] <- res$graph
res$graph <- tmp_graph

Amat <- res$graph

# force marker -> trait
Amat[1:num_phen, (num_phen + 1):num_var][Amat[1:num_phen, (num_phen + 1):num_var] != 0] <- 3
Amat[(num_phen + 1):num_var, 1:num_phen][Amat[(num_phen + 1):num_var, 1:num_phen] != 0] <- 2

writeMM(as(Amat, "sparseMatrix"), file = paste0(input_filestem, "_estimated_pag.mtx"))
print("Done")
