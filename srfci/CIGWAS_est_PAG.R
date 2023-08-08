library("pcalg");

source("/nfs/scistore13/robingrp/nmachnik/dev/ci-gwas/srfci/RFCI_functions_final.R")

myargs = commandArgs(trailingOnly=TRUE)
print(myargs)
input_filestem = myargs[1]
num_individuals = as.numeric(myargs[2])
output_file = myargs[3]
srfci_mode = myargs[4] # one of [mpu, mpd, std]
# sRFCI can be run in three different modes:
# mpu: limit unshielded triples to traits, leave marker-trait links undirected, run all subsequent steps as in original RFCI, orient marker-trait edges in the end as marker -> trait
# mpd: limit unshielded triples to traits, then direct marker-trait links as marker -> trait, then orient v-structures and run all subsequent steps
# std: run standard RFCI on the whole graph, without any modifications.

print(paste0("input_filestem: ", input_filestem))
print(paste0("num_individuals: ", num_individuals))
print(paste0("output_file: ", output_file))
print(paste0("mode: ", srfci_mode))

if (srfci_mode == "mpu") {
    unsh_triple_pheno_only = TRUE
    force_marker_to_trait_before_R1 = FALSE
    force_marker_to_trait_in_the_end = TRUE
} else if (srfci_mode == "mpd") {
    unsh_triple_pheno_only = TRUE
    force_marker_to_trait_before_R1 = TRUE
    force_marker_to_trait_in_the_end = FALSE
} else if (srfci_mode == "std") {
    unsh_triple_pheno_only = FALSE
    force_marker_to_trait_before_R1 = FALSE
    force_marker_to_trait_in_the_end = FALSE
} else {
    print("mode has to be one of [mpu, mpd, std]")
    exit()
}

form_sepset <- function(seps, nvar, max_level){
  res <- rep(list(rep(list(NULL), times=nvar)), nvar)
  lapply(1:nvar, function(i)
    lapply(1:nvar, function(j) 
      lapply(1:max_level, function(k) {
        val <- seps[ (((i - 1) * nvar + j) - 1) * max_level + k ]
        if (val != -1) {
          if (is.null(res[[i]][[j]])) {
            res[[i]][[j]] <<- c(val + 1)
          } else {
            res[[i]][[j]] <<- c(res[[i]][[j]], val + 1)
          }
        }
      })))
  
  return(res)
}

loadSparseSepSetIntoNestedList = function(filepath, nvar) {
  res <- rep(list(rep(list(NULL), times=nvar)), nvar)
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    vals = as.integer(strsplit(line, " +")[[1]])
    gvals <<- vals
    res[[vals[1]]][[vals[2]]] = tail(vals, -2)
  }
  close(con)
  return(res)
}

# number of variables, number of phenotypes, max_level
mdim_path <- paste0(input_filestem, ".mdim")
# separation set
sep_path <- paste0(input_filestem, ".ssm")

x <- read.csv(mdim_path, sep="\t", header=FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
max_level <- as.numeric(x[3])
sepset = loadSparseSepSetIntoNestedList(sep_path, num_var)

adjmat <- readMM(paste0(input_filestem, "_sam.mtx"))
cormat <- readMM(paste0(input_filestem, "_scm.mtx"))

suffStat <- list(C=cormat, n=num_individuals)
indepTest = gaussCItest
alpha=0.05
conservative = FALSE 
maj.rule = FALSE

print("Searching for unshielded triples")
if (unsh_triple_pheno_only) {
    u.t <- find.unsh.triple(adjmat[1:num_phen,1:num_phen], check = FALSE)
} else {
    u.t <- find.unsh.triple(adjmat, check = FALSE)
}

print("Orienting v-structures")
r.v. <- rfci.vStruc(suffStat=suffStat,
                    indepTest=gaussCItest,
                    alpha,
                    sepset,
                    adjmat, 
                    unshTripl=u.t$unshTripl,
                    unshVect=u.t$unshVect,
                    conservative=(conservative || maj.rule),
                    version.unf=c(1, 1),
                    maj.rule=maj.rule,
                    verbose=FALSE)
A <- r.v.$amat
sepset <- r.v.$sepset
if (force_marker_to_trait_before_R1) {
    A[1:num_phen,(num_phen+1):num_var][A[1:num_phen,(num_phen+1):num_var]==1]<-3
    A[(num_phen+1):num_var,1:num_phen][A[(num_phen+1):num_var,1:num_phen]==1]<-2
}

print("Applying R1-R10")
estimate_pag_with_traits_only = FALSE
if (estimate_pag_with_traits_only) {
    res <- udag2apag(A[1:num_phen,1:num_phen], suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, rules = rep(TRUE, 10), 
                 unfVect = r.v.$unfTripl, verbose = FALSE)
} else {
    res <- udag2apag(A, suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, rules = rep(TRUE, 10), 
                 unfVect = r.v.$unfTripl, verbose = FALSE)
}

Amat <- res$graph
MMfile3 <- file.path(output_file)

if (force_marker_to_trait_in_the_end) {
    Amat[1:num_phen,(num_phen+1):num_var][Amat[1:num_phen,(num_phen+1):num_var]==1]<-3
    Amat[(num_phen+1):num_var,1:num_phen][Amat[(num_phen+1):num_var,1:num_phen]==1]<-2
}

writeMM(Amat, file=MMfile3)
print("Done")
