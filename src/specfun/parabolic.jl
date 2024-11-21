# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Parabolic cylinder functions"""
#=

- ✅ pbvv
    - ✅ vvsa
    - ✅ vvla
- ✅ pbdv
    - ✅ dvsa
    - ✅ dvla
- ✅ pbwa
- cpbdn
    - ✅ cpdla
    - cpdsa
    - [gaih]
=#

"""
Compute parabolic cylinder functions W(a,±x)
and their derivatives

Input
- `a` --- Parameter  ( 0 ≤ |a| ≤ 5 )
- `x` --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )

Output
- `(w1f, w1d, w2f, w2d)`
    - `W1F` --- W(a,x)
    - `W1D` --- W'(a,x)
    - `W2F` --- W(a,-x)
    - `W2D` --- W'(a,-x)

Routine called:
- [`Specfun.cgama`](@ref) for computing complex gamma function
"""
function pbwa(a::Float64, x::Float64)
    EPS = 1e-15
    p0 = 0.59460355750136
    h = Vector{Float64}(undef, 100)
    d = Vector{Float64}(undef, 80)

    if a == 0.0
        g1 = 3.625609908222
        g2 = 1.225416702465
    else
        g1 = abs(cgama(0.25 + 0.5im * a, 1))
        g2 = abs(cgama(0.75 + 0.5im * a, 1))
    end

    f1 = sqrt(g1 / g2)
    f2 = sqrt(2.0 * g2 / g1)

    h0 = 1.0
    h1 = a
    h[1] = a
    for L1 in 2:100
        hl = a * h1 - 0.25 * (2 * L1 - 2.0) * (2 * L1 - 3.0) * h0
        h[L1] = hl
        h0, h1 = h1, hl
    end

    y1f = 1.0
    r = 1.0
    for k in 1:100
        r = 0.5 * r * x * x / (k * (2.0 * k - 1.0))
        r1 = h[k] * r
        y1f += r1
        if abs(r1) <= EPS * abs(y1f) && k > 30
            break
        end
    end

    y1d = a
    r = 1.0
    for k in 1:99
        r = 0.5 * r * x * x / (k * (2.0 * k + 1.0))
        r1 = h[k+1] * r
        y1d += r1
        if abs(r1) <= EPS * abs(y1d) && k > 30
            break
        end
    end
    y1d *= x

    d1 = 1.0
    d2 = a
    d[1] = 1.0
    d[2] = a
    for L2 in 3:80
        m = (2 * L2 - 1)
        dl = a * d2 - 0.25 * (m - 2.0) * (m - 3.0) * d1
        d[L2] = dl
        d1, d2 = d2, dl
    end

    y2f = 1.0
    r = 1.0
    for k in 1:79
        r = 0.5 * r * x * x / (k * (2.0 * k + 1.0))
        r1 = d[k+1] * r
        y2f += r1
        if abs(r1) <= EPS * abs(y2f) && k > 30
            break
        end
    end
    y2f *= x

    y2d = 1.0
    r = 1.0
    for k in 1:79
        r = 0.5 * r * x * x / (k * (2.0 * k - 1.0))
        r1 = d[k+1] * r
        y2d += r1
        if abs(r1) <= EPS * abs(y2d) && k > 30
            break
        end
    end

    w1f = p0 * (f1 * y1f - f2 * y2f)
    w2f = p0 * (f1 * y1f + f2 * y2f)
    w1d = p0 * (f1 * y1d - f2 * y2d)
    w2d = p0 * (f1 * y1d + f2 * y2d)
    return w1f, w1d, w2f, w2d
end


"""
Compute parabolic cylinder function Vv(x)
for large argument

Input
- `x`  --- Argument
- `va` --- Order

Output
- Vv(x)

Routines called:
- [`Specfun.dvla`](@ref) for computing Dv(x) for large |x|
- [`Specfun.gamma2`](@ref) for computing Г(x)
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
- `x`  --- Argument
- `va` --- Order

Output
- Dv(x)

Routines called:
- [`Specfun.vvla`](@ref) for computing Vv(x) for large |x|
- [`Specfun.gamma2`](@ref) for computing Г(x)
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
- `x`  --- Argument
- `va` --- Order

Output
- Dv(x)

Routine called
- [`Specfun.gamma2`](@ref) for computing Г(x)
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
- `x` --- Argument of Dv(x)
- `v` --- Order of Dv(x)

Output
- `DV(na)` --- Dn+v0(x)
- `DP(na)` --- Dn+v0'(x)
    ( na = |n|, v0 = v-n, |v0| < 1,
    n = 0,±1,±2,… )
- `(pdf, pdd)`
    - PDF --- Dv(x)
    - PDD --- Dv'(x)

Routines called:
- [`Specfun.dvsa`](@ref) for computing Dv(x) for small |x|
- [`Specfun.dvla`](@ref) for computing Dv(x) for large |x|
"""
function pbdv(dv::Vector{Float64}, dp::Vector{Float64}, x::Float64, v::Float64)
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
- `x`  --- Argument
- `va` --- Order

Output
- Vv(x)

Routine called
- [`Specfun.gamma2`](@ref) for computing Г(x)
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

"""
Compute parabolic cylinder functions Vv(x)
and their derivatives

Input
- `x` --- Argument of Vv(x)
- `v` --- Order of Vv(x)

Output
- `VV(na)` --- Vv(x)
- `VP(na)` --- Vv'(x)
    ( na = |n|, v = n+v0, |v0| < 1
     n = 0,±1,±2,… )
- `(pvf, pvd)`
    - `PVF` --- Vv(x)
    - `PVD` --- Vv'(x)

Routines called:
- [`Specfun.vvsa`](@ref) for computing Vv(x) for small |x|
- [`Specfun.vvla`](@ref) for computing Vv(x) for large |x|
"""
function pbvv(vv::Vector{Float64}, vp::Vector{Float64}, x::Float64, v::Float64)
    xa = abs(x)
    v += copysign(1.0, v)
    nv = trunc(Int, v)
    v0 = v - nv
    na = abs(nv)
    qe = exp(0.25 * x * x)
    q2p = sqrt(2.0 / pi)
    ja = na >= 1 ? 1 : 0

    # NOTE: When v==-0.0, access `vv[3]``
    arr_len = max(na+1, 3)
    # NOTE: in fortran, the index of `vv,vp` start from 0
    @assert length(vv) >= arr_len
    @assert length(vp) >= arr_len
    if v <= 0.0
        if v0 == 0.0
            if xa <= 7.5
                pv0 = vvsa(x, v0)
            else
                pv0 = vvla(x, v0)
            end
            f0 = q2p * qe
            f1 = x * f0
            vv[1] = pv0
            vv[2] = f0
            vv[3] = f1
        else
            for l in 0:ja
                v1 = v0 - l
                if xa <= 7.5
                    f1 = vvsa(x, v1)
                else
                    f1 = vvla(x, v1)
                end
                f0 = l == 0 ? f1 : f0
            end
            vv[1] = f0
            vv[2] = f1
        end

        kv = v0 == 0.0 ? 3 : 2
        for k in kv:na
            f = x * f1 + (k - v0 - 2.0) * f0
            vv[k+1] = f
            f0 = f1
            f1 = f
        end
    else
        if 0.0 ≤ x ≤ 7.5
            v2 = v < 1.0 ? (v + 1.0) : v
            f1 = vvsa(x, v2)
            v1 = v2 - 1.0
            kv = trunc(Int, v2)
            f0 = vvsa(x, v1)
            vv[kv+1] = f1
            vv[kv] = f0
            for k in (kv-2):-1:0
                f = x * f0 - (k + v0 + 2.0) * f1
                if k ≤ na
                    vv[k+1] = f
                end
                f1 = f0
                f0 = f
            end
        elseif x > 7.5
            pv0 = vvla(x, v0)
            m = 100 + abs(na)
            vv[2] = pv0
            f1 = 0.0
            f0 = 1.0e-40
            f = 0.0
            for k in m:-1:0
                f = x * f0 - (k + v0 + 2.0) * f1
                if k ≤ na
                    vv[k+1] = f
                end
                f1 = f0
                f0 = f
            end
            s0 = pv0 / f
            for k in 0:na
                vv[k+1] *= s0
            end
        else
            if xa ≤ 7.5
                f0 = vvsa(x, v0)
                v1 = v0 + 1.0
                f1 = vvsa(x, v1)
            else
                f0 = vvla(x, v0)
                v1 = v0 + 1.0
                f1 = vvla(x, v1)
            end
            
            vv[1] = f0
            vv[2] = f1
            for k in 2:na
                f = (x * f1 - f0) / (k + v0)
                vv[k+1] = f
                f0 = f1
                f1 = f
            end
        end
    end

    for k in 0:(na-1)
        v1 = v0 + k
        if v ≥ 0.0
            vp[k+1] =  0.5 * x * vv[k+1] - (v1 + 1.0) * vv[k+2]
        else
            vp[k+1] = -0.5 * x * vv[k+1] + vv[k+2]
        end
    end

    pvf = vv[na]
    pvd = vp[na]
    return pvf, pvd
end

"""
Compute complex parabolic cylinder function Dn(z)
for large argument

Input:
- z   --- Complex argument of Dn(z), Re(z) > 0, |z| >> n, z -> +Inf 
- n   --- Order of Dn(z) (n = 0, ±1, ±2,…)

Output:
- cdn --- Dn(z)
"""
function cpdla(n::Int, z::Complex{Float64})
    _EPS = 1e-12
    cb0 = z^n * exp(-0.25 * z*z)

    cr = 1.0 + 0im
    cdn = 1.0 + 0im
    for k in 1:16
        cr = -0.5 * cr * (2*k - n - 1) * (2*k - n - 2) / (k*z*z)
        cdn += cr
        if abs(cr) < abs(cdn) * _EPS
            break
        end
    end

    return cdn * cb0
end

"""
Compute complex parabolic cylinder function Dn(z)
for small argument

Input:
- z   --- Complex argument of Dn(z), Re(z) <= 3
- n   --- Order of Dn(z) (n = 0, -1, -2,…)

Output:
- cdn --- Dn(z)

Routine called:
- [`Specfun.gaih`](@ref) for computing Г(x), x = n/2 (n = 1, 2, ...)
"""
function cpdsa(n::Int, z::Complex{Float64})
    @assert real(z) <= 3
    @assert (-n) >= 0
    _EPS = 1.0e-15
    cdn = 0.0 + 0im

    ca0 = exp(-0.25 * z*z)
    va0 = 0.5 * (1.0 - n)
    if n == 0
        cdn = ca0
    elseif abs(z) == 0.0
        if (va0 <= 0.0) && (va0 == trunc(Int, va0))
            # TODO: unreachable, when n <= 0
            cdn = 0.0 + 0im
        else
            ga0 = gaih(va0)
            pd = sqrt(pi) / (2.0^(-0.5 * n) * ga0)
            cdn = complex(pd, 0.0)
        end
    else
        xn = -n
        g1 = gaih(xn)
        cb0 = 2.0^(-0.5 * n - 1.0) * ca0 / g1
        vt = -0.5 * n
        g0 = gaih(vt)

        cdn = complex(g0, 0.0)
        cr = 1.0 + 0im
        for m in 1:250
            vm = 0.5 * (m - n)
            gm = gaih(vm)
            cr = - cr * sqrt(2.0) * z / m
            cdw = gm * cr
            cdn = cdw + cdn
            if abs(cdw) < abs(cdn) * _EPS
                break
            end
        end
        cdn = cb0 * cdn
    end

    return cdn
end
