# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Legendre functions

sph_harm(m, n, theta, phi)
Compute spherical harmonics.
+ sph_harmonic_unsafe -> sph_harmonic -> cephes_poch


## specfun
lpmv(m, v, x)
    Associated Legendre function of integer order and real degree.
    + pmv_wrap -> specfun_lpmv

clpmn(m, n, z)
    Associated Legendre function of the first kind for complex arguments.
    + specfun_clpmn

lpn(n, z)
    Legendre function of the first kind.
    + specfun_lpn

lqn(n, z)
    Legendre function of the second kind.
    + specfun_lpn

lpmn(m, n, z)
    Sequence of associated Legendre functions of the first kind.
    + specfun_lpmn

lqmn(m, n, z)
    Sequence of associated Legendre functions of the second kind.
    + specfun_lqmn
"""
