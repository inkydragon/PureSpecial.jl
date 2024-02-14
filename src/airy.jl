# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Airy functions
airy(z)
    Airy functions and their derivatives.
    + airy_wrap -> cephes_airy
    + cairy_wrap -> amos_airy,amos_biry

airye(z)
    Exponentially scaled Airy functions and their derivatives.
    + cairy_wrap_e -> amos_airy,amos_biry
    + cairy_wrap_e_real -> amos_airy,amos_biry

ai_zeros(nt)
    Compute nt zeros and values of the Airy function Ai and its derivative.
    + _specfun.airyzo(nt, 1) -> specfun_airyzo

bi_zeros(nt)
    Compute nt zeros and values of the Airy function Bi and its derivative.
    + _specfun.airyzo(nt, 2) -> specfun_airyzo

itairy(x)
    Integrals of Airy functions
    + itairy_wrap -> specfun_itairy
"""
