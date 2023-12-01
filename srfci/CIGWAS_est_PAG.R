library("pcalg")
library("fastmatch")
library("Matrix")

source("/nfs/scistore17/robingrp/nmachnik/dev/ci-gwas/srfci/RFCI_functions_final.R")

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)
input_filestem <- myargs[1]
alpha <- as.numeric(myargs[2])
num_individuals <- as.numeric(myargs[3])
srfci_mode <- myargs[4] # one of [mpu, mpd, std]
# sRFCI can be run in three different modes:
# mpu: limit unshielded triples to traits, leave marker-trait links undirected, run all subsequent steps as in original RFCI, orient marker-trait edges in the end as marker -> trait
# mpd: limit unshielded triples to traits, then direct marker-trait links as marker -> trait, then orient v-structures and run all subsequent steps
# std: run standard RFCI on the whole graph, without any modifications.

print(paste0("input_filestem: ", input_filestem))
print(paste0("num_individuals: ", num_individuals))
print(paste0("alpha: ", alpha))
print(paste0("mode: ", srfci_mode))

if (srfci_mode == "mpu") {
    unsh_triple_pheno_only <- TRUE
    force_marker_to_trait_before_R1 <- FALSE
    force_marker_to_trait_in_the_end <- TRUE
    external_ambiguous_triples <- FALSE
} else if (srfci_mode == "mpd") {
    unsh_triple_pheno_only <- TRUE
    force_marker_to_trait_before_R1 <- TRUE
    force_marker_to_trait_in_the_end <- FALSE
    external_ambiguous_triples <- FALSE
} else if (srfci_mode == "std") {
    unsh_triple_pheno_only <- FALSE
    force_marker_to_trait_before_R1 <- FALSE
    force_marker_to_trait_in_the_end <- FALSE
    external_ambiguous_triples <- FALSE
} else if (srfci_mode == "cusk2") {
    unsh_triple_pheno_only <- FALSE
    force_marker_to_trait_before_R1 <- FALSE
    force_marker_to_trait_in_the_end <- TRUE
    external_ambiguous_triples <- TRUE
} else if (srfci_mode == "cusk2") {
} else {
    print("mode has to be one of [mpu, mpd, std, cusk2]")
    exit()
}

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

x <- read.csv(mdim_path, sep = "\t", header = FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
max_level <- as.numeric(x[3])
num_atr <- as.numeric(x[4])
sepset <- load_sparse_sepsets(sep_path, num_var)

# adjmat <- readMM(paste0(input_filestem, "_sam.mtx"))
# cormat <- readMM(paste0(input_filestem, "_scm.mtx"))
adjmat <- as.matrix(readMM(paste0(input_filestem, "_sam.mtx")))
cormat <- as.matrix(readMM(paste0(input_filestem, "_scm.mtx")))

atr_file <- file(atr_path, "rb")
atr <- readBin(atr_file, integer(), size = 4, n = num_atr * 3)
atrmat <- matrix(atr, ncol = 3, byrow = TRUE)

suffStat <- list(C = cormat, n = num_individuals)
indepTest <- gaussCItest
conservative <- FALSE
maj.rule <- FALSE

print("Searching for unshielded triples")
if (unsh_triple_pheno_only) {
    u.t <- find.unsh.triple(adjmat[1:num_phen, 1:num_phen], check = FALSE)
} else {
    u.t <- find.unsh.triple(adjmat, check = FALSE)
}

print("Orienting v-structures")
r.v. <- rfci.vStruc(
    suffStat = suffStat,
    indepTest = gaussCItest,
    alpha,
    sepset,
    adjmat,
    unshTripl = u.t$unshTripl,
    unshVect = u.t$unshVect,
    conservative = (conservative || maj.rule),
    version.unf = c(1, 1),
    maj.rule = maj.rule,
    verbose = FALSE
)
A <- r.v.$amat
sepset <- r.v.$sepset
if (force_marker_to_trait_before_R1) {
    A[1:num_phen, (num_phen + 1):num_var][A[1:num_phen, (num_phen + 1):num_var] != 0] <- 3
    A[(num_phen + 1):num_var, 1:num_phen][A[(num_phen + 1):num_var, 1:num_phen] != 0] <- 2
}

# optionally mark unshielded triples a-b-c as unfaithful / ambiguous,
# if b is in maximal sepset, but not in minimal sepset.
if (external_ambiguous_triples) {
    p <- num_var
    for (i in 1:nrow(atrmat)) {
        # increment because of one based indexing
        x <- atrmat[i, 1] + 1
        y <- atrmat[i, 2] + 1
        z <- atrmat[i, 3] + 1
        r.v.$unfTripl <- c(r.v.$unfTripl, triple2numb(p, x, y, z))
    }
}

writeMM(as(A, "sparseMatrix"), file = paste0(input_filestem, sprintf("_estimated_vStruc_%s.mtx", srfci_mode)))
print("Done")

print("Applying R1-R10")
estimate_pag_with_traits_only <- TRUE
if (estimate_pag_with_traits_only) {
    res <- udag2apag_ci_gwas(A[1:num_phen,1:num_phen], sepset, rules=rep(TRUE, 10), unfVect = r.v.$unfTripl)
    tmp_graph <- A
    tmp_graph[1:num_phen, 1:num_phen] <- res$graph
    res$graph <- tmp_graph
} else {
    res <- udag2apag(A,
        suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, rules = rep(TRUE, 10),
        unfVect = r.v.$unfTripl, verbose = FALSE
    )
}

Amat <- res$graph

if (force_marker_to_trait_in_the_end) {
    Amat[1:num_phen, (num_phen + 1):num_var][Amat[1:num_phen, (num_phen + 1):num_var] != 0] <- 3
    Amat[(num_phen + 1):num_var, 1:num_phen][Amat[(num_phen + 1):num_var, 1:num_phen] != 0] <- 2
}


# writeMM(Amat, file=paste0(input_filestem, sprintf("_estimated_pag_%s.mtx", srfci_mode)))
writeMM(as(Amat, "sparseMatrix"), file = paste0(input_filestem, sprintf("_estimated_pag_%s.mtx", srfci_mode)))
print("Done")
