#!/usr/bin/env Rscript
library("MendelianRandomization")
library("Matrix")

args = commandArgs(trailingOnly=TRUE)
cusk_output_path = args[1]
num_samples = as.numeric(args[2])
# expecting FALSE, TRUE input here
rm_non_adjacent = as.logical(args[3])
# expecting FALSE, TRUE input here
use_ld = as.logical(args[4])
output_file = args[5]

# TODO: it would be better to use merged_blocks_sam.mtx for the trait-trait adjacency
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

for (outcome_ix in 1:num_trait) {
    # we drop the rows which are parents of the outcome
    outcome_eff_rows = which(mxp_adj[, outcome_ix] == 1)
    by = B[-outcome_eff_rows, outcome_ix]
    byse = SE[-outcome_eff_rows, outcome_ix]
    if (rm_non_adjacent) {
        tested_traits = which(pxp_adj[, outcome_ix] == 1)
        bx = B[-outcome_eff_rows, tested_traits, drop=FALSE]
        bxse = B[-outcome_eff_rows, tested_traits, drop=FALSE]
    } else {
        tested_traits = seq(1, num_trait)[-outcome_ix]
        bx = B[-outcome_eff_rows, -outcome_ix, drop=FALSE]
        bxse = SE[-outcome_eff_rows, -outcome_ix, drop=FALSE]
    }
    # make sure that we have more ivs than exposures
    sufficient_ivs = dim(bx)[1] > dim(bx)[2]
    if ((length(tested_traits) > 0) && (sufficient_ivs)) {
        if (use_ld) {
            sub_ld_mat = full_ld_mat[-outcome_eff_rows, -outcome_eff_rows]
            input = MendelianRandomization::mr_mvinput(bx=bx, bxse=bxse, by=by, byse=byse, correlation=sub_ld_mat)
        } else {
            input = MendelianRandomization::mr_mvinput(bx=bx, bxse=bxse, by=by, byse=byse)
        }
        res = MendelianRandomization::mr_mvivw(input, robust=TRUE)
    }
    for (exposure_ix in 1:num_trait) {
        if (exposure_ix == outcome_ix) {
            next
        } else if ((exposure_ix %in% tested_traits) && (sufficient_ivs)) {
            if (exposure_ix > outcome_ix) {
                mvivw_exp_ix = exposure_ix - 1
            } else {
                mvivw_exp_ix = exposure_ix
            }
            sources = c(sources, exposure_ix)
            sinks = c(sinks, outcome_ix)
            effects = c(effects, res$Estimate[mvivw_exp_ix])
            ps = c(ps, res$Pvalue[mvivw_exp_ix])
        } else {
            sources = c(sources, exposure_ix)
            sinks = c(sinks, outcome_ix)
            effects = c(effects, 0.0)
            ps = c(ps, 1.0)
        }
    }
}

output = data.frame(
    source=sources,
    sink=sinks,
    effect=effects,
    p=ps,
    iv_selection=rep("cusk", length(sources))
)

write.table(output, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
