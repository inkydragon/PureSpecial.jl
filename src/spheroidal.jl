# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Spheroidal wave functions
pro_ang1(m, n, c, x[, out])
Prolate spheroidal angular function of the first kind and its derivative
+ prolate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa

pro_rad1(m, n, c, x[, out])
Prolate spheroidal radial function of the first kind and its derivative
+ prolate_radial1_nocv_wrap -> specfun_segv,specfun_rswfp

pro_rad2(m, n, c, x[, out])
Prolate spheroidal radial function of the second kind and its derivative
+ prolate_radial2_nocv_wrap -> specfun_segv,specfun_rswfp

obl_ang1(m, n, c, x[, out])
Oblate spheroidal angular function of the first kind and its derivative
+ oblate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa

obl_rad1(m, n, c, x[, out])
Oblate spheroidal radial function of the first kind and its derivative
+ oblate_radial1_nocv_wrap -> specfun_segv,specfun_rswfo

obl_rad2(m, n, c, x[, out])
Oblate spheroidal radial function of the second kind and its derivative.
+ oblate_radial2_nocv_wrap -> specfun_segv,specfun_rswfo

pro_cv(m, n, c[, out])
Characteristic value of prolate spheroidal function
+ prolate_segv_wrap -> specfun_segv

obl_cv(m, n, c[, out])
Characteristic value of oblate spheroidal function
+ oblate_segv_wrap -> specfun_segv


pro_cv_seq(m, n, c)
Characteristic values for prolate spheroidal wave functions.
+ _specfun.segv

obl_cv_seq(m, n, c)
Characteristic values for oblate spheroidal wave functions.
+ _specfun.segv


pro_ang1_cv(m, n, c, cv, x[, out])
Prolate spheroidal angular function pro_ang1 for precomputed characteristic value
+ prolate_aswfa_wrap -> specfun_aswfa

pro_rad1_cv(m, n, c, cv, x[, out])
Prolate spheroidal radial function pro_rad1 for precomputed characteristic value
+ prolate_radial1_wrap -> specfun_rswfp

pro_rad2_cv(m, n, c, cv, x[, out])
Prolate spheroidal radial function pro_rad2 for precomputed characteristic value
+ prolate_radial2_wrap -> specfun_rswfp

obl_ang1_cv(m, n, c, cv, x[, out])
Oblate spheroidal angular function obl_ang1 for precomputed characteristic value
+ oblate_aswfa_wrap -> specfun_aswfa

obl_rad1_cv(m, n, c, cv, x[, out])
Oblate spheroidal radial function obl_rad1 for precomputed characteristic value
+ oblate_radial1_wrap -> specfun_rswfo

obl_rad2_cv(m, n, c, cv, x[, out])
Oblate spheroidal radial function obl_rad2 for precomputed characteristic
+ oblate_radial2_wrap -> specfun_rswfo
"""
