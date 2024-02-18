# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Gamma and related functions

## cephes
gamma(z)
    gamma function.
    + cephes_Gamma,cgamma

gammaln(x)
    Logarithm of the absolute value of the gamma function.
    + cephes_lgam

loggamma(z)
    Principal branch of the logarithm of the gamma function.
    + loggamma_real,loggamma -> special::loggamma -> cephes::lgam

gammasgn(x)
    Sign of the gamma function.
    + cephes_gammasgn

gammainc(a, x)
    Regularized lower incomplete gamma function.
    + cephes_igam

gammaincinv(a, y)
    Inverse to the regularized lower incomplete gamma function.
    + cephes_igami

gammaincc(a, x)
    Regularized upper incomplete gamma function.
    + cephes_igamc

gammainccinv(a, y)
    Inverse of the regularized upper incomplete gamma function.
    + cephes_igamci

beta(a, b)
    Beta function.
    + cephes_beta

betaln(a, b)
    Natural logarithm of absolute value of beta function.
    + cephes_lbeta


## boost
betainc(a, b, x)
    Regularized incomplete beta function.
    + ibeta_float,ibeta_double @ boost

betaincc(a, b, x)
    Complement of the regularized incomplete beta function.
    + ibetac_float,ibetac_double @ boost

betaincinv(a, b, y)
    Inverse of the regularized incomplete beta function.
    + ibeta_inv_float,ibeta_inv_double @ boost

betainccinv(a, b, y)
    Inverse of the complemented regularized incomplete beta function.
    + ibetac_inv_float,ibetac_inv_double @ boost


psi(z)
    The digamma function.
    + digamma,cdigamma -> special::digamma

digamma(z)
    The digamma function.
    + -> psi
    
rgamma(z)
    Reciprocal of the gamma function.
    + cephes_rgamma
    + crgamma -> special::rgamma

polygamma(n, x)
    Polygamma functions.
    + gamma,zeta,psi

multigammaln(a, d)
    Returns the log of multivariate gamma, also sometimes called the generalized gamma.
    + loggam

poch(z, m)
    Pochhammer symbol.
    + cephes_poch
"""
