"""Mendelian Randomization implementations with cuda-skeleton input"""

import numpy as np
from scipy.stats import norm


def _se_of_corr(r: np.array, n: int):
    return (1 - r * r) / np.sqrt(n - 2)


def _get_neighbors(adj: np.array, node_ix: int, num_phen: int):
    nb = (adj[node_ix, :] > 0) & (adj[:, node_ix] > 0)
    return np.where(nb[:num_phen])[0]


def _get_parent_markers(adj: np.array, node_ix: int, num_phen: int):
    snp_ixs = np.where(adj[node_ix, :] != 0)[0]
    return snp_ixs[snp_ixs >= num_phen]


def _mvivw(
    adj: np.array,
    corr: np.array,
    num_samples: int,
    node_ix: int,
    num_trait: int,
    random_effects=True,
):
    """
    Ported from
    https://github.com/cran/MendelianRandomization/blob/master/R/mr_mvivw-methods.R
    """
    phen_nb = []
    ivs = []
    for p in _get_neighbors(adj, node_ix, num_trait):
        parents = _get_parent_markers(adj, p, num_trait)
        if len(parents) > 0:
            phen_nb.append(p)
            ivs.append(parents)

    ivs = np.unique(np.concatenate(ivs))
    ld = corr[np.ix_(ivs, ivs)]

    bx = corr[np.ix_(ivs, phen_nb)]
    by = corr[np.ix_(ivs, [node_ix])]
    se_by = _se_of_corr(r=by, n=num_samples)
    omega = np.outer(se_by, se_by) * ld
    omega_inv = np.linalg.pinv(omega)
    beta = np.linalg.pinv(bx.T @ omega_inv @ bx) @ bx.T @ omega_inv @ by
    # fixed effects
    beta_se = np.sqrt(np.diag(np.linalg.pinv(bx.T @ omega_inv @ bx)))
    if random_effects:
        rse = by - bx @ beta
        beta_se *= np.max(
            np.sqrt((rse.T @ omega_inv @ rse) / (len(ivs) - len(phen_nb)))
        )
    p_value = 2 * norm.cdf(-np.abs(beta.flatten() / beta_se))
    return beta.flatten(), p_value, phen_nb


def cig_w_mvivw(
    adj: np.array,
    corr: np.array,
    num_samples: int,
    num_trait: int,
    random_effects=True,
):
    """Run MV-IVW on a ci-gwas skeleton

    Args:
        adj (np.array): trait + marker adjacency matrix
        corr (np.array): trait + marker correlation matrix
        num_samples (int): number of samples (individuals)
        num_trait (int): number of traits
        random_effects (bool, optional): Calculate MVIVW p-values using a random effects model. Defaults to True.

    Returns:
        _type_: _description_
    """
    pvals = np.ones((num_trait, num_trait))
    betas = np.zeros_like(pvals)
    for i in range(num_trait):
        res = _mvivw(
            adj=adj,
            corr=corr,
            num_samples=num_samples,
            node_ix=i,
            num_trait=num_trait,
            random_effects=random_effects,
        )
        for e, p, nix in zip(*res):
            pvals[nix, i] = p
            betas[nix, i] = e
    return betas, pvals


# def get_neighbors_non_bidirected(pag: np.array, node_ix: int, num_phen: int):
#     adj_non_bidirected = ((pag[node_ix, :] != 2) | (pag[:, node_ix] != 2)) & (
#         (pag[node_ix, :] > 0) & (pag[:, node_ix] > 0)
#     )
#     return np.where(adj_non_bidirected[:num_phen])[0]

# def mvivw(
#     pag: np.array,
#     corr: np.array,
#     ld_matrix: np.array,
#     num_samples: int,
#     node_ix: int,
#     num_phen: int,
#     exclude_bidirected=True,
#     limit_effect_size=True,
#     random_effects=True,
# ):
#     if exclude_bidirected:
#         nb_fn = get_neighbors_non_bidirected
#     else:
#         nb_fn = get_neighbors
#     phen_nb = []
#     ivs = []
#     for p in nb_fn(pag, node_ix, num_phen):
#         parents = get_parent_markers(pag, p, num_phen)
#         if len(parents) > 0:
#             phen_nb.append(p)
#             ivs.append(parents)
#     if len(phen_nb) == 0:
#         return [0], [1], []

#     if ld_matrix is None:
#         ivs = np.unique(np.concatenate(ivs))
#         ld = corr[np.ix_(ivs, ivs)]
#     else:
#         ivs = np.unique(np.concatenate(ivs)) - num_phen
#         ld = ld_matrix[np.ix_(ivs, ivs)]

#     bx = corr[np.ix_(ivs, phen_nb)]
#     by = corr[np.ix_(ivs, [node_ix])]
#     se_by = se_of_corr(r=by, n=num_samples)
#     omega = np.outer(se_by, se_by) * ld
#     omega_inv = np.linalg.pinv(omega)
#     beta = np.linalg.pinv(bx.T @ omega_inv @ bx) @ bx.T @ omega_inv @ by
#     # fixed effects
#     beta_se = np.sqrt(np.diag(np.linalg.pinv(bx.T @ omega_inv @ bx)))
#     if random_effects:
#         rse = by - bx @ beta
#         beta_se *= np.max(
#             np.sqrt((rse.T @ omega_inv @ rse) / (len(ivs) - len(phen_nb)))
#         )
#     p_value = 2 * norm.cdf(-np.abs(beta.flatten() / beta_se))
#     sig_beta = beta[beta < 0.05]
#     if limit_effect_size and np.any(np.abs(sig_beta.flatten()) > 1):
#         return [0], [1], []
#     return beta.flatten(), p_value, phen_nb


# def cig_w_mvivw(
#     pag: np.array,
#     corr: np.array,
#     ld_matrix: np.array,
#     num_samples: int,
#     num_phen: int,
#     exclude_bidirected=True,
#     limit_effect_size=True,
#     random_effects=True,
# ):
#     pvals = np.zeros((num_phen, num_phen))
#     betas = np.zeros_like(pvals)
#     for i in range(num_phen):
#         res = mvivw(
#             pag=pag,
#             corr=corr,
#             ld_matrix=ld_matrix,
#             num_samples=num_samples,
#             node_ix=i,
#             num_phen=num_phen,
#             exclude_bidirected=exclude_bidirected,
#             limit_effect_size=limit_effect_size,
#             random_effects=random_effects,
#         )
#         for e, p, nix in zip(*res):
#             pvals[nix, i] = p
#             betas[nix, i] = e
#     return betas, pvals
