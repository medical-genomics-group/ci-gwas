

library("pcalg")
library("Matrix");

source("./RFunctions/DAVS_funs.R")
source("./RFunctions/rmdag6.R")
source("./RFunctions/functions_rfci_davs.R")


myargs = commandArgs(trailingOnly=TRUE)
id= as.numeric(myargs[1])




####################################

#for e1 
#MAX_LEVEL <- 5
MAX_LEVEL <- 6

form_sepset <- function(seps, nvar){
  res <- rep(list(rep(list(NULL), times=nvar)), nvar)
  
  lapply(1:nvar, function(i)
    lapply(1:nvar, function(j) 
      lapply(1:MAX_LEVEL, function(k) {
        val <- seps[ (((i - 1) * nvar + j) - 1) * MAX_LEVEL + k ]
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



#for alpha=e8-e2
stem = sprintf("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e%s/all_merged",id)
#stem = sprintf("/nfs/scistore13/robingrp/human_data/cigwas_estonia/bdpc_d1_l6_a1e%s_wo_bp/all_merged",id)
#stem2 = sprintf("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e%s/all_merged_replace",id)
#for alpha=e1
#stem = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l5_a1e1/all_merged"
#for alpha=0.05
#stem = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a05/all_merged"



mdim_path <- paste0(stem, ".mdim")
sep_path <- paste0(stem, ".ssm")

x <- read.csv(mdim_path, sep="\t", header=FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
num_var
num_phen
sepset = loadSparseSepSetIntoNestedList(sep_path, num_var)

#for alpha=e8-e2
adjmat <- readMM(sprintf("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e%s/all_merged_sam.mtx",id))
cormat <- readMM(sprintf("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e%s/all_merged_scm.mtx",id))
#adjmat <- readMM(sprintf("/nfs/scistore13/robingrp/human_data/cigwas_estonia/bdpc_d1_l6_a1e%s_wo_bp/all_merged_sam.mtx",id))
#cormat <- readMM(sprintf("/nfs/scistore13/robingrp/human_data/cigwas_estonia/bdpc_d1_l6_a1e%s_wo_bp/all_merged_scm.mtx",id))
#for alpha=e1

#adjmat <- readMM("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l5_a1e1/all_merged_sam.mtx")
#cormat <- readMM("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l5_a1e1/all_merged_scm.mtx")

#for alpha=0.05

#adjmat <- readMM("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a05/all_merged_sam.mtx")
#cormat <- readMM("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a05/all_merged_scm.mtx")


suffStat <- list(C=cormat ,n=458747)
n=458747
indepTest = gaussCItest
alpha=0.05
#p=dim(cormat)[1]

conservative = FALSE 
maj.rule = FALSE

#Forcing PAG for marker-->phen and find unshielded triple for phen
#adjmat[1:num_phen,num_phen+1:num_var][adjmat[1:num_phen,num_phen+1:num_var]==1]<-2

print("PAG Started")
u.t <- find.unsh.triple(adjmat[1:num_phen,1:num_phen], check = FALSE)

#print("PAG Started")
#u.t <- find.unsh.triple(adjmat, check = FALSE)
print("rfci.vStruc Started")
r.v. <- rfci.vStruc(suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, adjmat, 
                    unshTripl = u.t$unshTripl, unshVect = u.t$unshVect, conservative = (conservative || 
                     maj.rule), version.unf = c(1, 1), maj.rule = maj.rule, 
                   verbose = FALSE)
A <- r.v.$amat
sepset <- r.v.$sepset
A[1:num_phen,(num_phen+1):num_var][A[1:num_phen,(num_phen+1):num_var]==1]<-3
A[(num_phen+1):num_var,1:num_phen][A[(num_phen+1):num_var,1:num_phen]==1]<-2
A
print("udag2apag Started")
res <- udag2apag(A[1:num_phen,1:num_phen], suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, rules = rep(TRUE, 10), 
                 unfVect = r.v.$unfTripl, verbose = FALSE)

Amat <- res$graph

Amat 

#MMfile3 <- file.path(sprintf("/nfs/scistore13/robingrp/human_data/cigwas_estonia/bdpc_d1_l6_a1e%s_wo_bp",id), "all_merged_pag_mk3_ns_lut.mtx")
MMfile3 <- file.path(sprintf("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e%s",id),"pheno_pag_mk3_ns.mtx")
#MMfile3                       
#for alpha=e1
#MMfile3 <- file.path("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l5_a1e1", "all_merged_pag_sk_mk3.mtx")
#MMfile3 <- file.path("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l5_a1e1", "all_merged_pag_unsh_triple_mk3.mtx")

#for alpha=0.05
#MMfile3 <- file.path("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a05", "all_merged_pag_sk_mk3.mtx")
#MMfile3 <- file.path("/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a05", "all_merged_pag_mk3_ns_lut.mtx")

writeMM(Amat, file=MMfile3)
print("Amat finished")



#save(Amat,file = paste( "./CIGWAS_est_PAG/CIGWAS_est_PAG_lut_e",toString(id), ".RData", sep = "" ))
######################################
