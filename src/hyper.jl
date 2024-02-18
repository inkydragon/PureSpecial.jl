# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Hypergeometric functions

hyp2f1(a, b, c, z)
    Gauss hypergeometric function 2F1(a, b; c; z)
    + hyp2f1 -> cephes_hyp2f1 @cephes/hyp2f1.c
    + hyp2f1_complex @_hyp2f1.pxd

hyp1f1(a, b, x)
    Confluent hypergeometric function 1F1.
    + hyp1f1_double -> hyp1f1_wrap -> <hyp1f1_wrap> -> boost::math::hypergeometric_1F1
    +               -> hyp1f1_wrap -> specfun_chgm
    + chyp1f1_wrap -> specfun_cchg

hyperu(a, b, x)
    Confluent hypergeometric function U
    + hyperu -> poch -> cephes_poch
             -> hypU_wrap -> specfun_chgu

hyp0f1(v, z)
    Confluent hypergeometric limit function 0F1.
    + _hyp0f1_real -> _hyp0f1_asy@cython
    + _hyp0f1_cmplx -> cbesi_wrap @amos_wrappers.c
                    -> cbesj_wrap @amos_wrappers.c
"""
