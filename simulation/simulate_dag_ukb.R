

cd /nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/sim_1


library("Matrix")
library("BEDMatrix")

gen_rand_dag_obs_snp <- function(
    Tr,
    nL,
    deg,
    prob_pleio,
    lo_mp,
    hi_mp,
    lo_pp,
    hi_pp,
    cc,
    outdir,
    blockfile,
    plinkfile,
    plinkexe
) {
	num <- nL + Tr
	prob2 <- deg / (Tr)
    if(prob2 > 1) {
        prob2 = 1
    }
    d1 <- read.table(blockfile, header=F)
    chr <- d1[d1$V1 == cc,]
    bim <- read.table(paste0(plinkfile,".bim"), header=F)
    b2 <- bim[bim$V1 == cc,]
    fam <- read.table(paste0(plinkfile,".fam"), header=F)
    n = nrow(fam)

    # 1:20 blocks contain causal variants
    # 5 per trait, per block
    # create adjacency matrix
    CVs <- c()
    for(i in 1:20){
		x1 = chr$V2[i] + 1
		x2 = chr$V3[i] + 1
		d1 <- as.character(b2$V2[x1])
		d2 <- as.character(b2$V2[x2])
		tmp <- b2[which(b2$V2 == d1):which(b2$V2 == d2),]
		# select a spacing that isn't exactly the edge
		sp <- seq(5,dim(tmp)[1],round(dim(tmp)[1]/5,0)-5)
		tCV_1 <- sample(tmp$V2[sp[1]:sp[2]],num, replace=FALSE)
		tCV_2 <- sample(tmp$V2[c(sp[2]+1):sp[3]],num, replace=FALSE)
		tCV_3 <- sample(tmp$V2[c(sp[3]+1):sp[4]],num, replace=FALSE)
		tCV_4 <- sample(tmp$V2[c(sp[4]+1):sp[5]],num, replace=FALSE)
		tCV_5 <- sample(tmp$V2[c(sp[5]+1):sp[6]],num, replace=FALSE)
		tmp <- t(matrix(rep(diag(num),5),num,5*num))
		rownames(tmp) <- c(tCV_1,tCV_2,tCV_3,tCV_4,tCV_5)
		CVs <- rbind(CVs,tmp)
	}

	#add pleiotropy to adjacency matrix
	SNP <- nrow(CVs)
	for (i in 1:SNP) {
		iv <- which(CVs[i,] == 1)
		n_iv <- ncol(CVs) - length(iv)
		my_rbin <- rbinom(n_iv, 1, prob_pleio)
		which_element <- which(my_rbin != 0)
		if(length(which_element) != 0){
			CVs[i,][-iv] <- my_rbin
		}
	}

	#create adjacency for phenotypes
	Ys <- matrix(0, num, num)
	for (i in 1:num) {
        my_rbin <- rbinom(num-i, 1, prob2)
        which_element <- which(my_rbin != 0)
        if (length(which_element) != 0) {
            Ys[i, i + which_element] <- 1
        }
    }

    # bind together
    # create empty data frame for phenos
    AA <- rbind(CVs,Ys)
    x <- data.frame(matrix(vector(), n, SNP+num))

    #simulate marker-trait causal effects
    for (i in 1:SNP) {
        descendants = which(AA[i,] == 1)
        nb <- length(descendants)
        b <- runif(nb, lo_mp, hi_mp) * sign(rnorm(nb))
        AA[i, descendants] <- b
    }

    #simulate trait-trait causal effects
    for (i in (SNP+1):(nrow(AA))) {
        descendants = which(AA[i, ] == 1)
        nb <- length(descendants)
        b <- runif(nb, lo_pp, hi_pp) * sign(rnorm(nb))
        AA[i, descendants] <- b
    }

    #extract SNP IDs use plink to extract data
    SIDs <- data.frame(rownames(CVs))
    write.table(SIDs, paste0(outdir,"CV_snp_ids.txt"),row.names = FALSE,col.names=FALSE,quote=FALSE)
    system(paste0(plinkexe," --memory 40000 --threads 10 --bfile ",plinkfile," --extract ",outdir,"CV_snp_ids.txt --make-bed --out ",outdir,"CV_data"))

    #read in marker data
    G <- BEDMatrix(paste0(outdir,"CV_data.bed"), simple_names=TRUE)

    # re-order marker data to match adjacency matrix
    G <- G[,SIDs[,1]]

    # fill missing genetic data to mean and scale to unit variance
    for (i in 1:ncol(G)){
    	g_tmp <- G[,i]
    	g_tmp[is.na(g_tmp)] <- mean(g_tmp, na.rm=TRUE)
    	g_tmp <- scale(g_tmp)
    	x[,i] <- g_tmp
    }
 
    # create individual-level phenotypic data
    for (i in 2:num) {
        anc <- which(AA[, i] != 0)
        wa <- x[, anc, drop = FALSE]
        wa <- model.matrix(as.formula(paste("~0+", paste(names(wa), collapse="+"))), data.frame(wa))
        b <- AA[anc, i]
        g <- (wa %*% b)
        x[, (SNP + i)] <- g + rnorm(n, 0, sqrt(1 - var(g)))
    }
    xout <- x[,(SNP + i):num]
    return(list(x=xout, A=AA))
}

args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
id = args[2]

set.seed(id)

lo_mp <- 0.0001
hi_mp <- 0.01
lo_pp <- 0.001
hi_pp <- 0.2
nL = 2
Tr = 10
num <- nL + Tr
cc <- 1
prob_pleio <- 0.2
deg = 3
prob2 <- deg / (Tr)
#outdir = paste0("/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/sim_1/")
blockfile = paste0("/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000.blocks")
plinkfile = paste0("/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced")
plinkexe = paste0("/nfs/scistore17/robingrp/mrobinso/software/plink2")

dag <- gen_rand_dag_obs_snp(Tr, nL, deg, prob_pleio, 
	lo_mp, hi_mp, lo_pp, hi_pp, 
	cc, outdir, blockfile, plinkfile, plinkexe)

# output causal effects
writeMM(dag$A, file=file.path(paste0(outdir,"true_causaleffects_seed",id,"_pp",prob_pleio,".mtx")))

dag_data_mat <- data.matrix(dag$x)
corrs <- cor(dag_data_mat)
writeMM(as(corrs, "sparseMatrix"), file=file.path(outdir, paste("corr_n", toString(n), "_SNP_", toString(SNP),"_it_",toString(id) ,".mtx", sep = "" )))
