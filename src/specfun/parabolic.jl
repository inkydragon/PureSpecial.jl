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


"""
Compute parabolic cylinder function Vv(x)
for large argument

Input
x  --- Argument
va --- Order

Output
PV --- Vv(x)

Routines called:
(1) DVLA for computing Dv(x) for large |x|
(2) GAMMA2 for computing Г(x)
"""
function vvla(x::Float64, va::Float64)
    EPS = 1.0e-12

    r = 1.0
    pv = 1.0
    for k in 1:18
        # CoSF 13.4.5; DLMF 12.9.2:  Asymptotic Expansions
        r = 0.5 * r * (2.0 * k + va - 1.0) * (2.0 * k + va) / (k * x^2)
        pv += r
        if abs(r / pv) < EPS
            break
        end
    end

    qe = exp(0.25 * x^2)
    a0 = abs(x)^(-va - 1) * sqrt(2.0 / pi) * qe
    pv *= a0
    if x < 0.0
        # CoSF 13.3.11:
        x1 = -x
        pdl = dvla(x1, va)
        gl = gamma2(-va)
        dsl = sin(pi * va) * sin(pi * va)
        pv = dsl * gl / pi * pdl - cos(pi * va) * pv
    end

    return pv
end

"""
Compute parabolic cylinder functions Dv(x)
for large argument

Input
x  --- Argument
va --- Order

Output
PD --- Dv(x)

Routines called:
(1) VVLA for computing Vv(x) for large |x|
(2) GAMMA2 for computing Г(x)
"""
function dvla(x::Float64, va::Float64)
    EPS = 1.0e-12

    r = 1.0
    pd = 1.0
    for k in 1:16
        # CoSF 13.4.2; DLMF 12.9.1:  Asymptotic Expansions
        r = -0.5 * r * (2.0 * k - va - 1.0) * (2.0 * k - va - 2.0) / (k * x^2)
        pd += r
        if abs(r / pd) < EPS
            break
        end
    end

    ep = exp(-0.25 * x^2)
    a0 = abs(x)^va * ep
    pd *= a0
    if x < 0.0
        # CoSF 13.3.10:
        x1 = -x
        vl = vvla(x1, va)
        gl = gamma2(-va)
        pd = pi * vl / gl + cos(pi * va) * pd
    end

    return pd
end
