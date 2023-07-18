"""Synthetic generation of a skeleton + covariance matrix pair.

Examples
--------
>>> sample_size = 10000
>>> n = 100
>>> alpha = 10 ** -5
>>> c = rand_cov_mat(n, cov_factor=0.3)
>>> p = pcorr_mat(c)
>>> p = np.clip(p, -0.999, 0.999)
>>> idp = independent(p, alpha, sample_size, n - 2)
>>> plt.imshow(~idp)
"""

import numpy as np
from scipy.stats import random_correlation
from scipy.stats import norm

def rand_cov_mat(ndim: int, factor=1) -> np.ndarray:
    """Generate a random covariance matrix.

    Args:
        ndim (int): number of variables
        factor (int, optional): off-diagonal elements of matrix are multiplied with this value. Defaults to 1.

    Returns:
        np.ndarray: ndim x ndim covariance matrix.
    """
    rng = np.random.default_rng()
    vals = rng.random(ndim)
    eigs = vals / np.sum(vals) * ndim
    cm = random_correlation.rvs(eigs)
    return cm * factor * (np.eye(ndim) * 1 / factor  + np.ones_like(cm) - np.eye(ndim))

def pcorr_mat(covmat: np.ndarray) -> np.ndarray:
    """Compute partial correlation matrix.

    Partial correlations are computed as -(p_ij / sqrt(p_ii * p_jj)) where p is the inverse of covmat.

    Args:
        covmat (np.ndarray): covariance matrix

    Returns:
        np.ndarray: matrix of partial correlations
    """
    pmat = np.linalg.inv(covmat)
    ndim = pmat.shape[0]
    i, j = np.meshgrid(np.arange(ndim), np.arange(ndim))
    return pcorr(i, j, pmat)
    
def pcorr(i, j, pmat: np.ndarray):
    """Compute partial correlation of xi, xj | X \ {xi, xj}.

    Args:
        i (int or np.ndarray): _description_
        j (int or np.narray): _description_
        pmat (np.ndarray): _description_

    Returns:
        float or np.ndarray: partial correlation of xi, xj | X \ {xi, xj}
    """
    return -(pmat[i, j] / np.sqrt(pmat[i, i] * pmat[j, j]))

def fishers_z(pcorr):
    """Compute Fisher's z-transform of the partial correlation.

    Args:
        pcorr (float or np.ndarray): partial correlations

    Returns:
        float or np.ndarray: Fisher's z-transform of partial correlations
    """
    return np.abs(0.5 * np.log((1 + pcorr) / (1 - pcorr)))

def independent(pcorr, alpha, n, z) -> bool:
    """Test independence of two variables using their partial correlation.
    
    Tests H0: pcorr = 0 against two tail alternative Ha: pcorr != 0.
    Return true if the H0 is NOT rejected.

    Args:
        pcorr (float or np.ndarray): partial correlation
        alpha (float): significance level of H0
        n (int): sample size
        z (int): size of conditional set

    Returns:
        bool: True if the H0 is NOT rejected, i.e. if no evidence for dependence.
    """
    return (np.sqrt(n - z - 3) * fishers_z(pcorr)) <= norm.ppf(1 - alpha / 2)