library("pcalg");

source("/nfs/scistore13/robingrp/nmachnik/dev/ci-gwas/srfci/RFCI_functions_final.R")

myargs = commandArgs(trailingOnly=TRUE)
print(myargs)
input_filestem = myargs[1]
num_individuals = as.numeric(myargs[2])
output_file = myargs[3]

print(paste0("input_filestem: ", input_filestem))
print(paste0("num_individuals: ", num_individuals))
print(paste0("output_file: ", output_file))

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
unsh_triple_pheno_only = FALSE
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
                    verbose=TRUE)
A <- r.v.$amat
sepset <- r.v.$sepset
force_marker_to_trait = FALSE
if (force_marker_to_trait) {
    A[1:num_phen,(num_phen+1):num_var][A[1:num_phen,(num_phen+1):num_var]==1]<-3
    A[(num_phen+1):num_var,1:num_phen][A[(num_phen+1):num_var,1:num_phen]==1]<-2
}

print("Applying R5-R10")
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

writeMM(Amat, file=MMfile3)
print("Done")
