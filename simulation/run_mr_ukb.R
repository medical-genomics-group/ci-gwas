library(readr)
library(dplyr)
library(cause)
library(MRPRESSO)
library(MendelianRandomization)

myargs = commandArgs(trailingOnly=TRUE)
setseed = as.numeric(myargs[1])
sim = as.numeric(myargs[2])
outcome = as.numeric(myargs[3])
threshold = myargs[4]
eout = myargs[5]

if (is.na(setseed) || is.na(sim) || is.na(outcome) || is.na(threshold)) {
    stop("Received NA input argument, or argument is missing")
}

set.seed(setseed)
cause_res <- c()
mrpresso_res <- c()
mrpresso_mv_res <- c()
ivw_res <- c()
ivw_mv_res <- c()

indir <- paste("[DIR]")
infile <- paste0(indir,"sim_",sim,"/")
outf <- read.table(paste0(infile,"mxp_est.PHENO",outcome,".glm.linear"))

fseq <- seq(1,10,1)
fseq <- fseq[-outcome]

mf1 <- paste0(infile,"mxp_est.PHENO",fseq[1],".glm.linear")
mf2 <- paste0(infile,"mxp_est.PHENO",fseq[2],".glm.linear")
mf3 <- paste0(infile,"mxp_est.PHENO",fseq[3],".glm.linear")
mf4 <- paste0(infile,"mxp_est.PHENO",fseq[4],".glm.linear")
mf5 <- paste0(infile,"mxp_est.PHENO",fseq[5],".glm.linear")
mf6 <- paste0(infile,"mxp_est.PHENO",fseq[6],".glm.linear")
mf7 <- paste0(infile,"mxp_est.PHENO",fseq[7],".glm.linear")
mf8 <- paste0(infile,"mxp_est.PHENO",fseq[8],".glm.linear")
mf9 <- paste0(infile,"mxp_est.PHENO",fseq[9],".glm.linear")

mf_dat <- list(read.table(mf1),read.table(mf2),read.table(mf3),
                read.table(mf4),read.table(mf5),read.table(mf6),
                read.table(mf7),read.table(mf8),read.table(mf9))

snp1 <- read.table(paste0(mf1,".",threshold,".clumped.var"), header=TRUE)
snp2 <- read.table(paste0(mf2,".",threshold,".clumped.var"), header=TRUE)
snp3 <- read.table(paste0(mf3,".",threshold,".clumped.var"), header=TRUE)
snp4 <- read.table(paste0(mf4,".",threshold,".clumped.var"), header=TRUE)
snp5 <- read.table(paste0(mf5,".",threshold,".clumped.var"), header=TRUE)
snp6 <- read.table(paste0(mf6,".",threshold,".clumped.var"), header=TRUE)
snp7 <- read.table(paste0(mf7,".",threshold,".clumped.var"), header=TRUE)
snp8 <- read.table(paste0(mf8,".",threshold,".clumped.var"), header=TRUE)
snp9 <- read.table(paste0(mf9,".",threshold,".clumped.var"), header=TRUE)

top_vars <- list(snp1$SNP,snp2$SNP,snp3$SNP,snp4$SNP,snp5$SNP,
                snp6$SNP,snp7$SNP,snp8$SNP,snp9$SNP)

## Run CAUSE
for(i in 1:9){
    X <- cbind(mf_dat[[i]][,c(3,9,10,12)],outf[,c(9,10,12)])
    names(X) <- c("snp","beta_hat_1","seb1","p1",
                "beta_hat_2","seb2","p2")
    X <- new_cause_data(X)
    params <- est_cause_params(X, X$snp)
    res <- summary(cause(X=X, variants = top_vars[[i]], param_ests = params))
    v = c(paste0("y",outcome),paste0("y",fseq[i]), res$quants[[2]][1,1], res$p)
    cause_res = rbind(cause_res,v)
}
cause_res <- data.frame(cause_res)
cause_res <- cause_res[,c(2,1,4,3)]
names(cause_res) <- c("Exposure", "Outcome","p","est")

## Run MRPRESSO
for(i in 1:9){
    X <- cbind(outf[,c(3,9,10,12)],mf_dat[[i]][,c(9,10,12)])
    names(X) <- c("snp","beta_hat_1","seb1","p1",
                "beta_hat_2","seb2","p2")
    Xtop <- X[X$snp %in% top_vars[[i]],]
    mrp <- mr_presso(BetaOutcome = "beta_hat_1", 
                    BetaExposure = "beta_hat_2", 
                    SdOutcome = "seb1", 
                    SdExposure = "seb2", 
                    OUTLIERtest = TRUE, 
                    DISTORTIONtest = TRUE, 
                    data = Xtop, 
                    NbDistribution = 1000,  SignifThreshold = 0.1)
    if (is.na(mrp$'Main MR results'[2,"P-value"])) {
        v = c(paste0("y",outcome),paste0("y",fseq[i]), unlist(mrp$'Main MR results'[1,]))
        } else {
            v = c(paste0("y",outcome),paste0("y",fseq[i]), unlist(mrp$'Main MR results'[2,]))
            }
    mrpresso_res = rbind(mrpresso_res,v)
}
mrpresso_res <- data.frame(mrpresso_res)
mrpresso_res <- mrpresso_res[,c(2,1,8,5)]
names(mrpresso_res) <- c("Exposure", "Outcome","p","est")

## Run MV-MRPRESSO
X <- cbind(outf[,c(3,9,10,12)],
            mf_dat[[1]][,c(9,10,12)],
            mf_dat[[2]][,c(9,10,12)],
                mf_dat[[3]][,c(9,10,12)],
                mf_dat[[4]][,c(9,10,12)],
                mf_dat[[5]][,c(9,10,12)],
                mf_dat[[6]][,c(9,10,12)],
                mf_dat[[7]][,c(9,10,12)],
                mf_dat[[8]][,c(9,10,12)],
                mf_dat[[9]][,c(9,10,12)])
names(X) <- c("snp","beta_hat_1","seb1","p1",
                "beta_hat_2","seb2","p2",
                "beta_hat_3","seb3","p3",
                "beta_hat_4","seb4","p4",
                "beta_hat_5","seb5","p5",
                "beta_hat_6","seb6","p6",
                "beta_hat_7","seb7","p7",
                "beta_hat_8","seb8","p8",
                "beta_hat_9","seb9","p9",
                "beta_hat_10","seb10","p10")
ctvar <- unique(unlist(top_vars))
Xtop <- X[X$snp %in% ctvar,]
mrp <- mr_presso(BetaOutcome = "beta_hat_1", 
                BetaExposure = c("beta_hat_2", "beta_hat_3", "beta_hat_4",
                                 "beta_hat_5", "beta_hat_6", "beta_hat_7",
                                 "beta_hat_8", "beta_hat_9", "beta_hat_10"), 
                SdOutcome = "seb1", 
                SdExposure = c("seb2", "seb3", "seb4", "seb5",
                                "seb6", "seb7", "seb8", "seb9", "seb10"), 
                OUTLIERtest = TRUE, 
                DISTORTIONtest = TRUE, 
                data = Xtop, 
                NbDistribution = 1000,  
                SignifThreshold = 0.1)
if (is.na(mrp$'Main MR results'[mrp$'Main MR results'[,2] == "Outlier-corrected",][1,6])) {
    v = data.frame("v1" = paste0("y",outcome),
                    "v2" = paste0("y",fseq),
                    mrp$'Main MR results'[mrp$'Main MR results'[,2] == "Raw",])
    } else {
         v = data.frame("v1" = paste0("y",outcome),
                        "v2" = paste0("y",fseq),
                mrp$'Main MR results'[mrp$'Main MR results'[,2] == "Outlier-corrected",])
        }
mrpresso_mv_res = v
mrpresso_mv_res <- data.frame(mrpresso_mv_res)
mrpresso_mv_res <- mrpresso_mv_res[,c(2,1,8,5)]
names(mrpresso_mv_res) <- c("Exposure", "Outcome","p","est")


## Run IVW
for(i in 1:9){
X <- cbind(mf_dat[[i]][,c(3,9,10,12)],outf[,c(9,10,12)])
    names(X) <- c("snp","beta_hat_1","seb1","p1",
                "beta_hat_2","seb2","p2")
Xtop <- X[X$snp %in% top_vars[[i]],]
ivw <- mr_ivw(mr_input(bx = Xtop$beta_hat_1, bxse = Xtop$seb1, by = Xtop$beta_hat_2, byse = Xtop$seb2),
      robust = TRUE)
    v = data.frame("Exposure" = paste0("y",fseq[i]),
                    "Outcome" = paste0("y",outcome),
                    "p" = ivw$Pvalue,
                    "est" = ivw$Estimate,
                    "NumIVs" = ivw$SNPs)
ivw_res = rbind(ivw_res,v)
}


## Run MV-IVW
X <- cbind(outf[,c(3,9,10,12)],
            mf_dat[[1]][,c(9,10,12)],
            mf_dat[[2]][,c(9,10,12)],
                mf_dat[[3]][,c(9,10,12)],
                mf_dat[[4]][,c(9,10,12)],
                mf_dat[[5]][,c(9,10,12)],
                mf_dat[[6]][,c(9,10,12)],
                mf_dat[[7]][,c(9,10,12)],
                mf_dat[[8]][,c(9,10,12)],
                mf_dat[[9]][,c(9,10,12)])
names(X) <- c("snp","beta_hat_1","seb1","p1",
                "beta_hat_2","seb2","p2",
                "beta_hat_3","seb3","p3",
                "beta_hat_4","seb4","p4",
                "beta_hat_5","seb5","p5",
                "beta_hat_6","seb6","p6",
                "beta_hat_7","seb7","p7",
                "beta_hat_8","seb8","p8",
                "beta_hat_9","seb9","p9",
                "beta_hat_10","seb10","p10")
ctvar <- unique(unlist(top_vars))
Xtop <- X[X$snp %in% ctvar,]
ivw <- mr_mvivw(mr_mvinput(bx = cbind(Xtop$beta_hat_2, Xtop$beta_hat_3, Xtop$beta_hat_4,
                                    Xtop$beta_hat_5, Xtop$beta_hat_6, Xtop$beta_hat_7,
                                    Xtop$beta_hat_8, Xtop$beta_hat_9, Xtop$beta_hat_10),
                        bxse = cbind(Xtop$seb2, Xtop$seb3, Xtop$seb4, Xtop$seb5, Xtop$seb6, 
                            Xtop$seb7, Xtop$seb8, Xtop$seb9, Xtop$seb10), 
                        by = Xtop$beta_hat_1, byse = Xtop$seb1), robust = TRUE)
ivw_mv_res = data.frame("Exposure" = paste0("y",fseq),
                "Outcome" = paste0("y",outcome),    
                    "p" = ivw$Pvalue,
                    "est" = ivw$Estimate,
                    "NumIVs" = ivw$SNPs)

outdir <- paste0("[DIR]sim",sim,"_e",eout,"/")

write.csv(ivw_mv_res, paste0(outdir,"mr_mvivw_alpha",threshold,"_sim",sim,"_outcome_y",outcome,"_seed",setseed))
write.csv(ivw_res, paste0(outdir,"mr_ivw_alpha",threshold,"_sim",sim,"_outcome_y",outcome,"_seed",setseed))
write.csv(cause_res, paste0(outdir,"mr_cause_alpha",threshold,"_sim",sim,"_outcome_y",outcome,"_seed",setseed))
write.csv(mrpresso_res, paste0(outdir,"mr_presso_alpha",threshold,"_sim",sim,"_outcome_y",outcome,"_seed",setseed))
write.csv(mrpresso_mv_res, paste0(outdir,"mr_mvpresso_alpha",threshold,"_sim",sim,"_outcome_y",outcome,"_seed",setseed))
