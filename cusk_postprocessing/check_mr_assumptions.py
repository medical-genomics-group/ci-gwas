"""Fltering of exposures and instrument variables according to MR assumptions
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.io import mmread

def _fisher_z(v):
    return np.abs(0.5 * np.log(np.abs((1 + v) / (1 - v))))

def _alpha_thr(alpha: int, n: int, l: int):
    return norm.ppf(1 - (alpha / 2)) / np.sqrt(n - l - 3)

def _indep(
    x_ix: int,
    y_ix: int,
    s_ixs: list[int],
    corr: np.array,
    sample_size: int,
    alpha: float
):
    all_ixs = np.concatenate([[x_ix, y_ix], s_ixs]).astype(int)
    corr_sub = corr[np.ix_(all_ixs, all_ixs)]
    m = np.linalg.inv(corr_sub)
    pcorr = _fisher_z(-(m[0, 1] / np.sqrt(np.abs(m[0, 0] * m[1, 1]))))
    return pcorr < _alpha_thr(alpha, sample_size, len(s_ixs))

def _get_snp_parents(trait_ix, adj, num_traits):
    parents = np.where(adj[trait_ix, :])[0]
    return parents[parents > num_traits]

def _load_mdim(basepath: str):
    with open(basepath + ".mdim", "r") as fin:
        fields = fin.readline().strip().split("\t")
    return [int(f) for f in fields]


def get_iv_candidates(result_basename: str) -> pd.DataFrame:
    # get correlation files
    adj = mmread(f"{result_basename}_sam.mtx").toarray()
    _, num_traits, _ = _load_mdim(result_basename)
    trait_ixs = np.arange(0, num_traits)
    iv_candidates = {trait_ix: set(_get_snp_parents(trait_ix, adj, num_traits)) for trait_ix in trait_ixs}
    iv_snps = {(e, o): iv_candidates[e] for e in trait_ixs for o in trait_ixs if e != o}
    rows = []
    for (e_ix, o_ix), ivs in iv_snps.items():
        for snp_ix in ivs:
            # +1 for 1-based indexing in R
            rows.append({
                "Exposure": e_ix + 1,
                "Outcome": o_ix + 1,
                "IV": snp_ix + 1 - num_traits,
            })
    return pd.DataFrame(rows)


def check_ivs(
    result_basename: str,
    sample_size: int,
    accept_alpha: float,
    reject_alpha: float,
    relaxed_local_faithfulness=False,
    check_reverse_causality=False,
) -> pd.DataFrame:
    # get correlation files
    adj = mmread(f"{result_basename}_sam.mtx").toarray()
    pearson_corr = mmread(f"{result_basename}_scm.mtx").toarray()
    np.fill_diagonal(pearson_corr, 1)
    _, num_traits, _ = _load_mdim(result_basename)
    # num_snp = num_var - num_traits
    trait_ixs = np.arange(0, num_traits)
    # snp_ixs = np.arange(num_traits, num_var)
    iv_candidates = {trait_ix: set(_get_snp_parents(trait_ix, adj, num_traits)) for trait_ix in trait_ixs}

    # where do we have evidence for reverse causation?
    # outcome -> exposures
    rev_cause = {trait_ix: set() for trait_ix in trait_ixs}
    if check_reverse_causality:
        for outcome in trait_ixs:
            for exposure in trait_ixs:
                if exposure == outcome:
                    continue
                for snp in iv_candidates[outcome]:
                    marginal_dep = not _indep(snp, exposure, [], pearson_corr, sample_size, accept_alpha)
                    cond_indep = _indep(snp, exposure, [outcome], pearson_corr, sample_size, reject_alpha)
                    if marginal_dep and cond_indep:
                        rev_cause[outcome].add(exposure)

    trait_ixs_set = set(trait_ixs)
    valid_exposures = {trait_ix: trait_ixs_set - (rev_cause[trait_ix] | set([trait_ix])) for trait_ix in trait_ixs}

    # filter snps for valid ivs
    iv_snps = {(e, o): set() for e in trait_ixs for o in trait_ixs if e != o}
    for outcome in trait_ixs:
        for exposure in valid_exposures[outcome]:
            for snp in iv_candidates[exposure]:
                if relaxed_local_faithfulness:
                    # by assumption
                    marginal_dep = True
                else:
                    marginal_dep = not _indep(snp, outcome, [], pearson_corr, sample_size, accept_alpha)
                cond_indep = _indep(snp, outcome, list(valid_exposures[outcome]), pearson_corr, sample_size, reject_alpha)
                if cond_indep and marginal_dep:
                    iv_snps[(exposure, outcome)].add(snp)

    rows = []
    for (e_ix, o_ix), ivs in iv_snps.items():
        for snp_ix in ivs:
            # +1 for 1-based indexing in R
            rows.append({
                "Exposure": e_ix + 1,
                "Outcome": o_ix + 1,
                "IV": snp_ix + 1 - num_traits,
            })

    return pd.DataFrame(rows)
