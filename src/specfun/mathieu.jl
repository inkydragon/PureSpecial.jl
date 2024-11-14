# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Mathieu Functions"""
#=
mathieu_a(m, q)
mathieu_b(m, q)
+ specfun_cva2
    -> refine -> cvf
    -> cv0
    -> cvqm,cvql

mathieu_even_coef(m, q)
mathieu_odd_coef(m, q)
+ specfun_fcoef

mathieu_cem(m, q, x)
mathieu_sem(m, q, x)
+ specfun_mtu0
    -> cva2
    -> cva2

mathieu_modcem1(m, q, x)
mathieu_modcem2(m, q, x)
mathieu_modsem1(m, q, x)
mathieu_modsem2(m, q, x)
+ specfun_mtu12
    -> cva2,fcoef
    -> jynb -> jynbh
=#
