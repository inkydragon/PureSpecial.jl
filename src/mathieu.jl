# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Mathieu and related functions

mathieu_a(m, q)
    Characteristic value of even Mathieu functions
    + cem_cva_wrap -> specfun_cva2

mathieu_b(m, q)
    Characteristic value of odd Mathieu functions
    + sem_cva_wrap -> specfun_cva2

mathieu_even_coef(m, q)
    Fourier coefficients for even Mathieu and modified Mathieu functions.
    + _specfun.fcoef(kd=1,2) -> specfun_fcoef

mathieu_odd_coef(m, q)
    Fourier coefficients for even Mathieu and modified Mathieu functions.
    + _specfun.fcoef(kd=3,4) -> specfun_fcoef

mathieu_cem(m, q, x)
    Even Mathieu function and its derivative
    + cem_wrap -> specfun_mtu0

mathieu_sem(m, q, x)
    Odd Mathieu function and its derivative
    + sem_wrap -> specfun_mtu0

mathieu_modcem1(m, q, x)
    Even modified Mathieu function of the first kind and its derivative
    + mcm1_wrap -> specfun_mtu12

mathieu_modcem2(m, q, x)
    Even modified Mathieu function of the second kind and its derivative
    + mcm2_wrap -> specfun_mtu12

mathieu_modsem1(m, q, x)
    Odd modified Mathieu function of the first kind and its derivative
    + msm1_wrap -> specfun_mtu12

mathieu_modsem2(m, q, x)
    Odd modified Mathieu function of the second kind and its derivative
    + msm2_wrap -> specfun_mtu12
"""
