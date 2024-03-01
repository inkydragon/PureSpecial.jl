# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Parabolic Cylinder functions

pbwa(a, x)
    Parabolic cylinder function W
    pbwa_wrap -> specfun_pbwa -> specfun_cgama

pbvv(v, x)
    Parabolic cylinder function V
    pbvv_wrap -> specfun_pbvv
                    + specfun_vvsa -> gamma2
                    + specfun_vvla
    
pbdv(v, x)
    Parabolic cylinder function D
    pbdv_wrap -> specfun_pbdv
                    + specfun_dvsa -> specfun_vvla -> gamma2
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

"""
Compute parabolic cylinder function Dv(x)
for small argument

Input
x  --- Argument
va --- Order

Output
PD --- Dv(x)

Routine called
GAMMA2 for computing Г(x)
"""
function dvsa(x::Float64, va::Float64)
    EPS = 1.0e-12
    sq2 = sqrt(2.0)

    ep = exp(-0.25 * x^2)
    va0 = 0.5 * (1.0 - va)
    pd = NaN
    if va == 0.0
        pd = ep
    else
        if x == 0.0
            if va0 <= 0.0 && isinteger(va0)
                pd = 0.0
            else
                ga0 = gamma2(va0)
                pd = sqrt(pi) / (2.0^(-0.5 * va) * ga0)
            end
        else
            g1 = gamma2(-va)
            a0 = 2.0^(-0.5 * va - 1.0) * ep / g1
            vt = -0.5 * va
            g0 = gamma2(vt)
            pd = g0
            r = 1.0
            for m in 1:250
                vm = 0.5 * (m - va)
                gm = gamma2(vm)
                r = -r * sq2 * x / m
                r1 = gm * r
                pd += r1
                if abs(r1) < abs(pd) * EPS
                    break
                end
            end
            pd *= a0
        end
    end

    return pd
end

"""
Compute parabolic cylinder functions Dv(x)
and their derivatives

Input
x --- Argument of Dv(x)
v --- Order of Dv(x)

Output
DV(na) --- Dn+v0(x)
DP(na) --- Dn+v0'(x)
    ( na = |n|, v0 = v-n, |v0| < 1,
    n = 0,±1,±2,… )
PDF --- Dv(x)
PDD --- Dv'(x)

Routines called:
(1) DVSA for computing Dv(x) for small |x|
(2) DVLA for computing Dv(x) for large |x|
"""
function pbdv!(dv::Vector{Float64}, dp::Vector{Float64}, x::Float64, v::Float64)
    xa = abs(x)
    v += copysign(1.0, v)
    nv = trunc(Int, v)
    v0 = v - nv
    na = abs(nv)
    ep = exp(-0.25 * x * x)
    ja = na >= 1 ? 1 : 0

    # NOTE: in fortran, the index of `dv,dp` start from 0
    @assert length(dv) >= (na+1)
    @assert length(dp) >= (na+1)
    if v >= 0.0
        if v0 == 0.0
            pd0 = ep
            pd1 = x * ep
        else
            for l in 0:ja
                v1 = v0 + l
                if xa <= 5.8
                    pd1 = dvsa(x, v1)
                else
                    pd1 = dvla(x, v1)
                end
                pd0 = l == 0 ? pd1 : pd0
            end
        end
        dv[1] = pd0
        dv[2] = pd1
        for k in 2:na
            pdf = x * pd1 - (k + v0 - 1.0) * pd0
            dv[k+1] = pdf
            pd0 = pd1
            pd1 = pdf
        end
    else
        if x <= 0.0
            if xa <= 5.8
                pd0 = dvsa(x, v0)
                v1 = v0 - 1.0
                pd1 = dvsa(x, v1)
            else
                pd0 = dvla(x, v0)
                v1 = v0 - 1.0
                pd1 = dvla(x, v1)
            end
            dv[1] = pd0
            dv[2] = pd1
            for k in 2:na
                pd = (-x * pd1 + pd0) / (k - 1.0 - v0)
                dv[k+1] = pd
                pd0 = pd1
                pd1 = pd
            end
        elseif x <= 2.0
            v2 = nv + v0
            if nv == 0
                v2 -= 1.0
            end
            nk = trunc(Int, -v2)
            f1 = dvsa(x, v2)
            v1 = v2 + 1.0
            f0 = dvsa(x, v1)
            dv[nk+1] = f1
            dv[nk] = f0
            for k in (nk-2):-1:0
                f = x * f0 + (k - v0 + 1.0) * f1
                dv[k+1] = f
                f1 = f0
                f0 = f
            end
        else
            if xa <= 5.8
                pd0 = dvsa(x, v0)
            else
                pd0 = dvla(x, v0)
            end
            dv[1] = pd0
            m = 100 + na
            f1 = 0.0
            f0 = 1e-30
            f = 0.0
            for k in m:-1:0
                f = x * f0 + (k - v0 + 1.0) * f1
                if k <= na
                    dv[k+1] = f
                end
                f1 = f0
                f0 = f
            end
            s0 = pd0 / f
            for k in 0:na
                dv[k+1] *= s0
            end
        end
    end

    for k in 0:(na-1)
        v1 = abs(v0) + k
        if v >= 0.0
            dp[k+1] =  0.5 * x * dv[k+1] - dv[k+2]
        else
            dp[k+1] = -0.5 * x * dv[k+1] - v1 * dv[k+2]
        end
    end

    pdf = dv[na]
    pdd = dp[na]
    return pdf, pdd
end


"""
Compute parabolic cylinder function Vv(x)
for small argument

Input
x  --- Argument
va --- Order

Output
PV --- Vv(x)

Routine called
GAMMA2 for computing Г(x)
"""
function vvsa(x::Float64, va::Float64)
    EPS = 1.0e-12
    ep = exp(-0.25 * x * x)
    va0 = 1.0 + 0.5 * va

    if x == 0.0
        if (va0 <= 0.0 && isinteger(va0)) || va == 0.0
            pv = 0.0
        else
            vb0 = -0.5 * va
            sv0 = sin(va0 * pi)
            ga0 = gamma2(va0)
            pv = 2.0^vb0 * sv0 / ga0
        end
    else
        sq2 = sqrt(2.0)
        a0 = 2.0^(-0.5 * va) * ep / (2.0 * pi)
        sv = sin(-(va + 0.5) * pi)
        v1 = -0.5 * va
        g1 = gamma2(v1)
        pv = (sv + 1.0) * g1
        r = 1.0
        fac = 1.0

        for m in 1:250
            vm = 0.5 * (m - va)
            gm = gamma2(vm)
            r = r * sq2 * x / m
            fac = -fac
            gw = fac * sv + 1.0
            r1 = gw * r * gm
            pv += r1
            if (abs(r1 / pv) < EPS) && (gw != 0.0)
                break
            end
        end
        pv *= a0
    end

    return pv
end
