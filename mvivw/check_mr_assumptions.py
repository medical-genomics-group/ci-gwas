## load mxp, pxp, corrs and ess, put into one matrix

"""
All that we want to do here is to check MR assumptions:

A) Avoid reverse causation
    Filter out exposures for which we find a marker that satisfies
        1) marker dep exposure | {}
        2) marker indep exposure | outcome

B) Valid IVs
    Filter out markers which don't satisy
        1) marker dep outcome | {}
        2) marker indep outcome | exposures
"""

from scipy.stats import norm
from scipy.io import mmread

def fisher_z(v):
    return np.abs(0.5 * np.log(np.abs((1 + v) / (1 - v))))

def alpha_thr(alpha: int, n: int, l: int):
    return norm.ppf(1 - (alpha / 2)) / np.sqrt(n - l - 3)

def indep(x_ix, y_ix, S_ixs, corr, ess, alpha):
    all_ixs = np.concatenate([[x_ix, y_ix], S_ixs]).astype(int)
    ess_sub = ess[np.ix_(all_ixs, all_ixs)]
    corr_sub = corr[np.ix_(all_ixs, all_ixs)]
    ess_mean = np.nanmean(ess_sub[np.triu_indices(len(all_ixs), 1)])
    m = np.linalg.inv(corr_sub)
    pcorr = fisher_z(-(m[0, 1] / np.sqrt(np.abs(m[0, 0] * m[1, 1]))))
    return pcorr < alpha_thr(alpha, ess_mean, len(S_ixs))

def get_snp_parents(trait_ix, adj, num_traits):
    parents = np.where(adj[trait_ix, :])[0]
    return parents[parents > num_traits]

wdir = "/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/mr_local_faithfulness"
pxp_c_file = f"{wdir}/pxp_corrs.tsv"
mxp_c_file = f"{wdir}/mxp_corrs.tsv"
pxp_ess_file = f"{wdir}/pxp_ess.tsv"
mxp_ess_file = f"{wdir}/mxp_ess.tsv"

mxp_ess = pd.read_csv(mxp_ess_file, sep="\t")
pxp_ess = pd.read_csv(pxp_ess_file, sep="\t")
mxp_c = pd.read_csv(mxp_c_file, sep="\t")
pxp_c = pd.read_csv(pxp_c_file, sep="\t")

n_traits = pxp_ess.shape[0]
n_snps = mxp_ess.shape[0]
n_var = n_traits + n_snps

trait_ixs = np.arange(0, n_traits)
snp_ixs = np.arange(n_traits, n_var)

corr = np.zeros(shape=(n_var, n_var))
ess = np.zeros(shape=(n_var, n_var))

corr[np.ix_(trait_ixs, trait_ixs)] = pxp_c
corr[np.ix_(trait_ixs, snp_ixs)] = mxp_c.T
corr[np.ix_(snp_ixs, trait_ixs)] = mxp_c

diag = np.arange(n_var)
corr[diag, diag] = 1

ess[np.ix_(trait_ixs, trait_ixs)] = pxp_ess
ess[np.ix_(trait_ixs, snp_ixs)] = mxp_ess.T
ess[np.ix_(snp_ixs, trait_ixs)] = mxp_ess

trait_names = mxp_ess.columns

adj = mmread("/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/ukb_8m/cuskss/all_alpha_e4/cuskss_merged_sam.mtx").toarray()
pearson_corr = mmread("/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/ukb_8m/cuskss/all_alpha_e4/cuskss_merged_scm.mtx").toarray()
full_ss_matrix = np.ones_like(ess) * 458747

binary_traits = ["CAD", "AT", "T2D", "SMK", "ST"]
binary_trait_ixs = [np.where(trait_names == trait_name)[0][0] for trait_name in binary_traits]
iv_candidates = {trait_ix: set(get_snp_parents(trait_ix, adj, n_traits)) for trait_ix in trait_ixs}

alpha = 0.0001

# we re-evaluate the rparents for binary exposures,
# because we found that there may be a high false positive rate in x -> y_cont -> y_bin situations when CI(x, y_bin | y_cont) is tested.

binary_traits = ["CAD", "AT", "T2D", "SMK", "ST"]
binary_trait_ixs = [np.where(trait_names == trait_name)[0][0] for trait_name in binary_traits]
iv_candidates = {trait_ix: get_snp_parents(trait_ix, adj, n_traits) for trait_ix in trait_ixs}

for outcome in trait_ixs:
    for exposure in binary_trait_ixs:
        if exposure == outcome:
            continue
        residual_ivs = set()
        for snp in iv_candidates[exposure]:
            if not indep(snp, exposure, [outcome], corr, ess, alpha):
                residual_ivs.add(snp)
        iv_candidates[exposure] = residual_ivs

# where do we have evidence for reverse causation?
# outcome -> exposures
rev_cause = {trait_ix: set() for trait_ix in trait_ixs}
for outcome in trait_ixs:
    for exposure in trait_ixs:
        if exposure == outcome:
            continue
        for snp in iv_candidates[outcome]:
            if exposure in binary_trait_ixs:
                marginal_dep = not indep(snp, exposure, [], corr, ess, alpha)
                cond_indep = indep(snp, exposure, [outcome], corr, ess, alpha)
            else:
                marginal_dep = not indep(snp, exposure, [], pearson_corr, full_ss_matrix, alpha)
                cond_indep = indep(snp, exposure, [outcome], pearson_corr, full_ss_matrix, alpha)
            if marginal_dep and cond_indep:
                rev_cause[outcome].add(exposure)

trait_ixs_set = set(trait_ixs)
valid_exposures = {trait_ix: trait_ixs_set - (rev_cause[trait_ix] | set([trait_ix])) for trait_ix in trait_ixs}

# filter snps for valid ivs
iv_snps = {(e, o): set() for e in trait_ixs for o in trait_ixs if e != o}
for exposure in trait_ixs:
    for outcome in trait_ixs:
        for snp in iv_candidates[exposure]:
            if outcome in binary_trait_ixs:
                cond_indep = indep(snp, outcome, list(valid_exposures[outcome]), corr, ess, alpha)
            else:
                cond_indep = indep(snp, outcome, list(valid_exposures[outcome]), pearson_corr, full_ss_matrix, alpha)
            if cond_indep:
                iv_snps[(exposure, outcome)].add(snp)