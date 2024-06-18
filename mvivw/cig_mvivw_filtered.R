#!/usr/bin/env Rscript
library("MendelianRandomization")
library("Matrix")

args = commandArgs(trailingOnly=TRUE)
cusk_output_path = args[1]
num_samples = as.numeric(args[2])
iv_df_path = args[3]
output_file = args[4]

iv_df = read.csv(iv_df_path)

corr_path = sprintf("%s/cuskss_merged_scm.mtx", cusk_output_path)
adj_path = sprintf("%s/cuskss_merged_sam.mtx", cusk_output_path)
mdim_path = sprintf("%s/cuskss_merged.mdim", cusk_output_path)

mdim = read.table(mdim_path, sep="\t")
num_trait = mdim[1, 2]
num_var = mdim[1, 1]
num_snp = num_var - num_trait
corrs = Matrix::readMM(corr_path)
adj = Matrix::readMM(adj_path)

full_ld_mat = data.matrix(corrs[(num_trait + 1):num_var, (num_trait + 1):num_var])
pxp_adj = adj[1:num_trait, 1:num_trait]
mxp_adj = data.matrix(t(adj[1:num_trait, (num_trait+1):num_var]))
B = data.matrix(t(corrs[1:num_trait, (num_trait+1):num_var]))
SE = (1 - B * B) / sqrt(num_samples - 2)

est_eff = matrix(0, num_trait, num_trait)
est_p = matrix(1.0, num_trait, num_trait)

# TODO: refactor, this is terrible
sources = c()
sinks = c()
effects = c()
ps = c()
sk_adj = c()

for (outcome_ix in 1:num_trait) {
    ivs = unique(iv_df[iv_df["Outcome"] == outcome_ix, ][, "IV"])
    exposures = unique(iv_df[iv_df["Outcome"] == outcome_ix, ][, "Exposure"])

    bx = B[ivs, exposures, drop=FALSE]
    bxse = SE[ivs, exposures, drop=FALSE]
    by = B[ivs, outcome_ix]
    byse = SE[ivs, outcome_ix]
    
    # make sure that we have more ivs than exposures
    sufficient_ivs = dim(bx)[1] > dim(bx)[2]
    if ((length(exposures) > 0) && (sufficient_ivs)) {
        input = MendelianRandomization::mr_mvinput(bx=bx, bxse=bxse, by=by, byse=byse)
        # res = MendelianRandomization::mr_mvivw(input, robust=TRUE)
        res = tryCatch(
            MendelianRandomization::mr_mvivw(input, robust=TRUE),
            warning = function(w) {
                print("Got a warning!")
                print(w)
                print(outcome_ix)
                res = MendelianRandomization::mr_mvivw(input, robust=TRUE)
                print(res)
                return(res)
            }
        )
    }
    for (exposure_ix in 1:num_trait) {
        if (exposure_ix == outcome_ix) {
            next
        } else if ((exposure_ix %in% exposures) && (sufficient_ivs)) {
            if (exposure_ix > outcome_ix) {
                mvivw_exp_ix = exposure_ix - 1
            } else {
                mvivw_exp_ix = exposure_ix
            }
            sources = c(sources, exposure_ix)
            sinks = c(sinks, outcome_ix)
            effects = c(effects, res$Estimate[mvivw_exp_ix])
            ps = c(ps, res$Pvalue[mvivw_exp_ix])
            sk_adj = c(sk_adj, pxp_adj[exposure_ix, outcome_ix] == 1)
        } else {
            sources = c(sources, exposure_ix)
            sinks = c(sinks, outcome_ix)
            effects = c(effects, 0.0)
            ps = c(ps, 1.0)
            sk_adj = c(sk_adj, pxp_adj[exposure_ix, outcome_ix] == 1)
        }
    }
}

output = data.frame(
    source=sources,
    sink=sinks,
    effect=effects,
    p=ps,
    sk_adj=sk_adj
)

write.table(output, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
