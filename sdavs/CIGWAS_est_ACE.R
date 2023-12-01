library("Matrix")
library("pcalg")
library("fastmatch")

source("/nfs/scistore17/robingrp/nmachnik/dev/ci-gwas/sdavs/DAVS_functions_final.R")

myargs <- commandArgs(trailingOnly = TRUE)
exposure <- as.numeric(myargs[1])
outcome <- as.numeric(myargs[2])
num_individuals <- as.numeric(myargs[3])
alpha <- as.numeric(myargs[4])
pag_path <- myargs[5]
skeleton_results_filestem <- myargs[6]
output_file <- myargs[7]

if (is.na(exposure) || is.na(outcome) || is.na(num_individuals) || is.na(pag_path) || is.na(skeleton_results_filestem) || is.na(output_file) || is.na(alpha)) {
    stop("Received NA input argument, or argument is missing")
}

max_k <- 3
max_depth <- NULL

mdim_path <- paste0(skeleton_results_filestem, ".mdim")
x <- read.csv(mdim_path, sep = "\t", header = FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
cormat <- readMM(paste0(skeleton_results_filestem, "_scm.mtx"))
pag <- readMM(sprintf(pag_path))
pag <- matrix(pag, nrow(pag), ncol(pag), dimnames = NULL)
cr.dat <- as(cormat, "matrix")

res_davs <- array(0, dim = c(num_var))
if (exposure != outcome) {
    wpa <- searchAM(pag, exposure, type = "pa")
    wsp <- searchAM(pag, exposure, type = "sp")
    ypa <- searchAM(pag, outcome, type = "pa")
    ysp <- searchAM(pag, outcome, type = "sp")
    ww <- union(wpa, wsp)
    yy <- union(ypa, ysp)
    QQ <- setdiff(ww, yy)
    for (q in QQ) {
        print(sprintf("DAVS at exposure: %s, outcome: %s, q: %s", exposure, outcome, q))
        res_davs[q] <- Davs.con.causaleffect_cor_new(cr.dat,
            exposure,
            outcome,
            q,
            pag,
            num_individuals,
            alpha = alpha,
            models = "linearreg",
            max_k = max_k,
            max_depth = max_depth,
            force_directed = FALSE
        )
    }
}

res_davs[res_davs == 0] <- NA
ace <- mean(res_davs, na.rm = TRUE)
write.csv(ace, file = file.path(output_file), row.names = FALSE)

warning()
