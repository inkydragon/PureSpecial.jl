# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Parabolic Cylinder functions

pbwa(a, x)
    Parabolic cylinder function W
    pbwa_wrap -> specfun_pbwa -> specfun_cgama

pbvv(v, x)
    Parabolic cylinder function V
    pbvv_wrap -> specfun_pbvv
                    + specfun_vvsa -> specfun_gamma2
                    + specfun_vvla
    
pbdv(v, x)
    Parabolic cylinder function D
    pbdv_wrap -> specfun_pbdv
                    + specfun_dvsa -> specfun_vvla -> specfun_gamma2
                    + specfun_dvla

"""
