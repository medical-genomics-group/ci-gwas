#!/usr/bin/env Rscript
library("MendelianRandomization")
library("Matrix")

args = commandArgs(trailingOnly=TRUE)
cusk_output_stem = args[1]
num_samples = as.numeric(args[2])
# expecting FALSE, TRUE input here
rm_non_adjacent = as.logical(args[3])
# expecting FALSE, TRUE input here
use_ld = as.logical(args[4])
rm_counterfactual = as.logical(args[5])
fixed_links_path = args[6]
output_file = args[7]

corr_path = sprintf("%s_scm.mtx", cusk_output_stem)
adj_path = sprintf("%s_sam.mtx", cusk_output_stem)
mdim_path = sprintf("%s.mdim", cusk_output_stem)

mdim = read.table(mdim_path, sep="\t")
num_trait = mdim[1, 2]
num_var = mdim[1, 1]
num_snp = num_var - num_trait
corrs = Matrix::readMM(corr_path)
adj = Matrix::readMM(adj_path)

if (rm_counterfactual) {
    fixed_links_file <- file(fixed_links_path, "rb")
    fixed_links <- readBin(fixed_links_file, integer(), size = 4, n = num_trait * num_trait)
    fixed_links <- matrix(fixed_links, ncol = num_trait, byrow = TRUE)
}

full_ld_mat = data.matrix(corrs[(num_trait + 1):num_var, (num_trait + 1):num_var])
pxp_adj = adj[1:num_trait, 1:num_trait]
mxp_adj = data.matrix(t(adj[1:num_trait, (num_trait+1):num_var]))
B = data.matrix(t(corrs[1:num_trait, (num_trait+1):num_var]))
SE = (1 - B * B) / sqrt(num_samples - 2)

est_eff = matrix(0, num_trait, num_trait)
est_p = matrix(1.0, num_trait, num_trait)

# TODO: refactor
sources = c()
sinks = c()
effects = c()
ps = c()
sk_adj = c()
num_snp = c()

for (outcome_ix in 1:num_trait) {
    # we drop the rows which are parents of the outcome
    outcome_eff_rows = which(mxp_adj[, outcome_ix] == 1)
    rix = seq(1, dim(B)[1])
    rix_keep = rix[! rix %in% outcome_eff_rows]
    if (rm_non_adjacent) {
        tested_traits = which(pxp_adj[, outcome_ix] == 1)
        # rm rows which are parents of the removed traits
        rm_rows = outcome_eff_rows
        for (rm_trait_ix in which(pxp_adj[, outcome_ix] != 1)) {
            rm_rows = c(rm_rows, which(mxp_adj[, rm_trait_ix] == 1))
        }
        bx = B[-rm_rows, tested_traits, drop=FALSE]
        bxse = B[-rm_rows, tested_traits, drop=FALSE]
        by = B[-rm_rows, outcome_ix]
        byse = SE[-rm_rows, outcome_ix]
    } else if (rm_counterfactual) {
        # rm traits which are definitely downstream of outcome
        rm_rows = which(fixed_links[outcome_ix, ] == 1)
        rm_rows = c(rm_rows, outcome_ix)
        rm_rows = c(rm_rows, outcome_eff_rows)
        tested_traits = which(fixed_links[outcome_ix, ] != 1)
        tested_traits = tested_traits[-outcome_ix]
        bx = B[-rm_rows, tested_traits, drop=FALSE]
        bxse = B[-rm_rows, tested_traits, drop=FALSE]
        by = B[-rm_rows, outcome_ix]
        byse = SE[-rm_rows, outcome_ix]
    } else {
        tested_traits = seq(1, num_trait)[-outcome_ix]
        bx = B[rix_keep, -outcome_ix, drop=FALSE]
        bxse = SE[rix_keep, -outcome_ix, drop=FALSE]
        by = B[rix_keep, outcome_ix]
        byse = SE[rix_keep, outcome_ix]
    }
    # make sure that we have more ivs than exposures
    sufficient_ivs = dim(bx)[1] > dim(bx)[2]
    if ((length(tested_traits) > 0) && (sufficient_ivs)) {
        if (use_ld) {
            sub_ld_mat = full_ld_mat[rix_keep, rix_keep]
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
            sk_adj = c(sk_adj, pxp_adj[exposure_ix, outcome_ix] == 1)
            num_snps = c(num_snp, dim(bx)[1])
        } else {
            sources = c(sources, exposure_ix)
            sinks = c(sinks, outcome_ix)
            effects = c(effects, 0.0)
            ps = c(ps, 1.0)
            sk_adj = c(sk_adj, pxp_adj[exposure_ix, outcome_ix] == 1)
            num_snps = c(num_snp, dim(bx)[1])
        }
    }
}

output = data.frame(
    source=sources,
    sink=sinks,
    effect=effects,
    p=ps,
    sk_adj=sk_adj,
    num_snps=num_snps
)

write.table(output, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
