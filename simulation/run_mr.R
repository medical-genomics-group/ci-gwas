# NM: do we really need all these packages?

library("igraph");library("bnlearn")
library("MRPRESSO");library("mixtools")
library("optparse");library("limma")
library("MendelianRandomization")
library("Matching");library("rmatio")
library("Zelig");library("Metrics")
library("CBPS");library("pcalg")
library("MatchIt");library("cobalt")
library("fastmatch");library("stdReg")
library("Matrix");library("readr")
library("dplyr");library("cause")
library("parallel");library("ParallelPC")
library("ff", quietly = TRUE)
library("ivreg");library("mclust")
library("tictoc")

run_lm <- function(
    x,
    y,
    z,
    df
) {
    if(is.null(z)){
        df = data.frame(x=df[,x],y=df[,y])
    } else{
        # this never happens, as we always specify z = NULL?
        df = data.frame(x=df[,x],y=df[,y],df[,z])
    }
    model = lm(x~.,data=df)
    coefs = summary(model)$coefficients
    return(coefs[2,])
}

run_cause_on_tr <- function(
    tr1,
    phenos,
    G_it,
    GWAS_effects,
    GWAS_ses
) {
    cause_res = c()
    for(tr2 in phenos) {
        if (tr1 == tr2){
            next
        }
        if (is.null(G_it)) {
            ivs = rownames(GWAS_effects)
        } else {
            ivs = rownames(G_it)[G_it[,tr1]]
        }
    
        X = cbind(
            GWAS_effects[ivs, tr1], GWAS_ses[ivs, tr1],
            GWAS_effects[ivs, tr2], GWAS_ses[ivs, tr2]
        )
    
        colnames(X) = c("beta_hat_1","seb1", "beta_hat_2","seb2")
        X = data.frame(X)
        X$snp = ivs
        X = new_cause_data(X)
        m_p = est_cause_params(X,variants = ivs)
        m = cause(X,m_p)
        # these p values are causal model better than null and causal better than shared.
        # we don't use them though; is the p value in the summary the right one?
        p1 = pnorm(m$elpd[2,5])
        p2 = pnorm(m$elpd[3,5]) 
        m_s = summary(m)
        # NM: @MR, if you could double check this, most importantly that m_s$p is the p-value
        # NM: for the causal model being better than the shared model..
        # NM: the tutorial is here: https://jean997.github.io/cause/ldl_cad.html
        v = c(tr1, tr2, m_s$quants[[2]][1,1], m_s$p)
        cause_res = rbind(cause_res, v)
  }
  
  return(cause_res)
}


args = commandArgs(trailingOnly=TRUE)
indir = args[1]
id = as.numeric(args[2])
n=16000
SNP=1600
dag_data = readRDS(file.path(indir, paste("dag_data_n",  toString(n), "_SNP_", toString(SNP),"_it_", toString(id),".rds", sep = "")))

all_mr_res = list()
df = dat
ivs = colnames(df)[grepl("X", colnames(df))]
phenos = colnames(df)[grepl("Y", colnames(df))]
num_ivs = length(ivs)
p = length(phenos)
# Get all IV-phenotype associations
GWAS_Ps = matrix(1, num_ivs, p, dimnames = list(ivs,phenos))
GWAS_effects = matrix(0, num_ivs, p, dimnames = list(ivs,phenos))
GWAS_ses = matrix(0, num_ivs, p, dimnames = list(ivs,phenos))
GWAS_Zs = matrix(0, num_ivs, p, dimnames = list(ivs,phenos))

for (pheno in phenos) {
    print(pheno)
    gwas_res = sapply(ivs, run_lm, x=pheno, z=NULL, df=df)
    GWAS_Ps[,pheno] = gwas_res[4,]
    GWAS_effects[,pheno] = gwas_res[1,]
    GWAS_ses[,pheno] = gwas_res[2,]
    GWAS_Zs[,pheno] = gwas_res[3,]
}

print(alpha_e)
G_it = GWAS_Ps < alpha_e
res_path="./mr_res_git_thr"

# NM: this is needed only for tpr / fpr for the 
# NM: marker-trait links, but we don't do that anymore
# MMfile_git <- file.path(res_path, paste("git_n", toString(n), "_SNP_", toString(SNP),  "_alpha_", toString(alpha_e), "_it_", toString(id), ".csv", sep = ""))
# print(MMfile_git)
# print(G_it)
# write.csv(G_it, file=MMfile_git,row.names = FALSE)

# NM: why is pleio size 100 a standard MR analysis?
# Run MR
# pleio size is set to 100 - no filtering of variants (a standard MR analysis)

egger_res <- run_pairwise_mr_analyses(
    G_it,
    GWAS_effects,
    GWAS_ses,
    pleio_size=100,
    pruned_lists=NULL,
    func=mr_egger,
    robust=T)

write.csv(
    egger_res,
    file=file.path(res_path, paste("mr_egger_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = "")),
    row.names=FALSE)

ivw_res <- run_pairwise_mr_analyses(
    G_it,
    GWAS_effects,
    GWAS_ses,
    pleio_size=100,
    pruned_lists=NULL,
    func=mr_ivw,
    robust=T) 

write.csv(
    ivw_res,
    file=file.path(res_path, paste("mr_ivw_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = "")),
    row.names=FALSE)

mrpresso_res = c()
# I don't know why these try blocks are here
try({
    for (tr1 in phenos) {
        currivs = rownames(GWAS_Ps)[G_it[,tr1]]
        for (tr2 in phenos) {
            if (tr1 == tr2) {
                next
            }
            X = data.frame(
                    E1b=GWAS_effects[currivs,tr1],
                    O1b=GWAS_effects[currivs,tr2],
                    E1sd=GWAS_ses[currivs,tr1],
                    O1sd=GWAS_ses[currivs,tr2])
            try({
                res = mr_presso(
                    BetaOutcome="O1b",
                    BetaExposure="E1b", 
                    SdOutcome="O1sd",
                    SdExposure="E1sd",
                    data=X,
                    OUTLIERtest=T,
                    DISTORTIONtest=T,
                    NbDistribution=1000,
                    SignifThreshold=0.1)
                if (is.na(res$`Main MR results`[2,"P-value"])) {
                    v = c(tr1, tr2, unlist(res$`Main MR results`[1,]))
                } else {
                    v = c(tr1, tr2, unlist(res$`Main MR results`[2,]))
                }
                v["GlobalTestP"] = res$`MR-PRESSO results`$`Global Test`$Pvalue
                mrpresso_res = rbind(mrpresso_res,v)
            })
        }
    }
    if (!is.null(dim(mrpresso_res))) {
        mrpresso_res = as.data.frame(mrpresso_res)
        for (j in 3:ncol(mrpresso_res)) {
            mrpresso_res[[j]] = as.numeric(as.character(mrpresso_res[[j]]))
        }
    }
})

write.csv(
    mrpresso_res
    file=file.path(res_path, paste("mr_mrpresso_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = ""))),
    row.names=FALSE)

cause_res = mclapply(
    phenos,
    run_cause_on_tr,
    phenos=phenos,
    G_it=G_it,
    GWAS_effects=GWAS_effects,
    GWAS_ses=GWAS_ses,
    mc.cores=4)

cause_res_all = c()
for (tr in cause_res) {
    if (length(tr) == 1) {
        next
    }
    cause_res_all = rbind(cause_res_all, tr)
}

cause_res_all_df = data.frame(cause_res_all, stringsAsFactors=F)

for(j in 3:ncol(cause_res_all_df)) {
    cause_res_all_df[[j]] = as.numeric(as.character(cause_res_all_df[[j]]))
}
 
write.csv(
    cause_res_all_df,
    file=file.path(res_path, paste("mr_cause_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = "")),
    row.names=FALSE)