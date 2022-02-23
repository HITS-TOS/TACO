# Multivariate effective sample size
# ==================================
#
# Reference https://research-engine.appspot.com/5898535211892736/bibliography/170001
#
# Define a *multivariate effective sample size* (mESS) as:
#
# mESS = n (det(Lambda)|det(Sigma))^{1/p}
#
# where n is the number of samples, Lambda is the covariance matrix of the MCMC chain,
# Sigma is the asymptotic covariance matrix in the *central limit theorem* (CLT)
# and p is the number of quantities estimated. The asymptotic matrix Sigma can be estimated
# by *batch means*. A stopping criteria can be defined when
#
# mESS >=  (2^{2/p} pi chi^2_{1-\alpha,p} ) / ([p Gamma(p/2)]^{2/p} epsilon^2)
#
# with Gamma being the gamma function, alpha the confidence estimate, chi^2 $ the chi square
# quantile function and epsilon is the tolerance accepted as the Monte-Carlo error.

import math
import scipy
import scipy.stats
import numpy as np


def batch_means(chain, batch_size):
    """
    Chain is assumed to be a numpy array of shape [n, p] with n being
    the total number of iterations and p the number of parameters.
    """
    n = chain.shape[0]
    num_batches = int(np.floor(n / batch_size))
    chain = chain[(n - (num_batches * batch_size)):, :]  # Reshape the chain so we have complete batches
    overall_mean = np.mean(chain, axis=0)
    _batch_means = np.zeros((num_batches, chain.shape[1]))
    for i in range(num_batches):
        batch_chain = chain[(i * batch_size):((i + 1) * batch_size), :]
        _batch_means[i, :] = np.mean(batch_chain, axis=0)
    return _batch_means - overall_mean


def multiESS(chain, nu=0.5):
    """
    Number of multivariate effective samples in chain with a batch means estimation
    of the posterior covariance structure with batches of size n^nu (where n is the
    number of iterations in the chain).
    """
    n = chain.shape[0]
    p = chain.shape[1]

    batch_size = int(np.floor(n ** nu))
    num_batches = int(np.floor(n / batch_size))
    chain = chain[(n - (num_batches * batch_size)):, :]  # Reshape the chain so we have complete batches
    n = chain.shape[0]
    lambda_mat = np.cov(chain, rowvar=False)
    sigma_mat = np.cov(batch_means(chain, batch_size), rowvar=False)
    lambda_det_p = np.linalg.det(lambda_mat) ** (1.0 / p)
    sigma_det_p = np.linalg.det(sigma_mat) ** (1.0 / p)
    mess = n * (lambda_det_p / sigma_det_p) / batch_size
    if np.isfinite(mess):
        return mess
    else:
        return 0


def minESS(p, alpha=0.05, eps=0.05):
    """
    Mimimum number of multivariate effective samples (integer) in order
    to have an estimate within an confidence interval 1-alpha with an error, 
    due to the Monte-Carlo estimation, of eps.
    """
    crit = scipy.stats.chi2.ppf(1 - alpha, p)
    foo = 2.0 / p
    logminESS = foo * math.log(2) + math.log(math.pi) - foo * math.log(p) - foo * scipy.special.gammaln(
        p / 2) - 2 * math.log(eps) + math.log(crit)
    return int(round(math.exp(logminESS)))
