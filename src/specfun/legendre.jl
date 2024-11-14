# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Legendre Functions

Note:
    - m isa Integer; m >= 0
    - n isa Integer; n >= 0

Legendre polynomial
- lpn:      Legendre polynomial, Pn(x), Pn'(x), x isa Real
- clpn:     Legendre polynomial, Pn(z), Pn'(z), z isa Complex, 
- lpni:     Legendre polynomial and its integral, Pn(x), Pn'(x), Integral of Pn(t) from 0 to x
- legzo:    Zeros of the Legendre polynomial Pn(x) in the interval [-1,1]

Legendre Function
- lqna:     Legendre function, Qn(x), Qn'(x), |x| <= 1
- lqnb:     Legendre function, Qn(x), Qn'(x), |x| > 1
- clqn:     Legendre function, Qn(z), Qn'(z)

Associated Legendre Function
- lpmn:     Associated legendre function, Pmn(x), Pmn'(x)
- clpmn:    Associated legendre function, Pmn(z), Pmn'(z)
- ✅ lpmns:    Associated legendre function for a given order, Pmn(x), Pmn'(x) 
- lqmn:     Associated legendre function, Qmn(x), Qmn'(x)
- clqmn:    Associated legendre function, Qmn(z), Qmn'(z)
- ✅ lqmns:    Associated legendre function for a given order, Qmn(x), Qmn'(x)
- ✅ lpmv:     Associated legendre function, Pmv(x), |x| <= 1, v isa Real, v >= 0, 
"""

"""
Compute associated Legendre functions Pmn(x)
and Pmn'(x) for a given order

Input:
- m --- Order of Pmn(x),  m = 0,1,2,...,n
- n --- Degree of Pmn(x), n = 0,1,2,...,N
- x --- Argument of Pmn(x)

Output:
- PM(n) --- Pmn(x)
- PD(n) --- Pmn'(x)
"""
function lpmns(m::Int, n::Int, x::Float64, pm::Vector{Float64}, pd::Vector{Float64})
    # NOTE: f77:  DIMENSION PM(0:N),PD(0:N)
    @assert m >= 0
    @assert n >= 0
    @assert length(pm) >= (n+1)
    @assert length(pm) >= 2
    @assert length(pm) >= (m+2)
    @assert length(pd) >= (n+1)

    for k in 0:n
        pm[k + 1] = 0.0
        pd[k + 1] = 0.0
    end

    if abs(x) == 1.0
        for k in 0:n
            if m == 0
                pm[k + 1] = 1.0
                pd[k + 1] = 0.5 * k * (k + 1.0)
                if x < 0.0
                    pm[k + 1] *= (k % 2 == 0 ? 1 : -1)
                    pd[k + 1] *= ((k + 1) % 2 == 0 ? 1 : -1)
                end
            elseif m == 1
                pd[k + 1] = 1e300
            elseif m == 2
                pd[k + 1] = -0.25 * (k + 2.0) * (k + 1.0) * k * (k - 1.0)
                if x < 0.0
                    pd[k + 1] *= ((k + 1) % 2 == 0 ? 1 : -1)
                end
            end
        end

        return
    end

    x0 = abs(1.0 - x * x)
    pm0 = 1.0
    pmk = pm0
    for k in 1:m
        pmk = (2.0 * k - 1.0) * sqrt(x0) * pm0
        pm0 = pmk
    end
    pm1 = (2.0 * m + 1.0) * x * pm0
    pm[m + 1] = pmk
    pm[m + 2] = pm1
    for k in (m+2):n
        pm2 = ((2.0 * k - 1.0) * x * pm1 - (k + m - 1.0) * pmk) / (k - m)
        pm[k + 1] = pm2
        pmk = pm1
        pm1 = pm2
    end

    pd[1] = ((1.0 - m) * pm[2] - x * pm[1]) / (x * x - 1.0)
    for k in 1:n
        pd[k + 1] = (k * x * pm[k + 1] - (k + m) * pm[k]) / (x * x - 1.0)
    end
    coef = ifelse(m % 2 == 0, 1, -1)
    for k in 1:n
        pm[k + 1] *= coef
        pd[k + 1] *= coef
    end
end

"""
Compute associated Legendre functions Qmn(x)
and Qmn'(x) for a given order

Input:
- m --- Order of Qmn(x),  m = 0,1,2,...
- n --- Degree of Qmn(x), n = 1,2,...
- x --- Argument of Qmn(x)

Output:
- QM(n) --- Qmn(x)
- QD(n) --- Qmn'(x)
"""
function lqmns(m::Int, n::Int, x::Float64, qm::Vector{Float64}, qd::Vector{Float64})
    # NOTE: f77:  DIMENSION QM(0:N),QD(0:N)
    @assert m >= 0
    @assert n >= 1
    @assert length(qm) >= (n+1)
    @assert length(qm) >= 2
    @assert length(qd) >= (n+1)

    val = ifelse(abs(x) == 1.0, 1e300, 0.0)
    for k in 0:n
        qm[k+1] = val
        qd[k+1] = val
    end
    if abs(x) == 1.0
        return
    end

    ls = ifelse(abs(x) > 1.0, -1, 1)
    xq = sqrt(ls * (1.0 - x * x))
    q0 = 0.5 * log(abs((x + 1.0) / (x - 1.0)))
    q00 = q0
    q10 = -1.0 / xq
    q01 = x * q0 - 1.0
    q11 = -ls * xq * (q0 + x / (1.0 - x * x))
    qf0 = q00
    qf1 = q10
    qm0 = 0.0
    qm1 = 0.0
    for k in 2:m
        qm0 = -2.0 * (k - 1.0) / xq * x * qf1 - ls * (k - 1.0) * (2.0 - k) * qf0
        qf0 = qf1
        qf1 = qm0
    end

    if m == 0
        qm0 = q00
    elseif m == 1
        qm0 = q10
    end
    qm[1] = qm0
    if abs(x) < 1.0001
        if (m == 0) && (n > 0)
            qf0 = q00
            qf1 = q01
            for k in 2:n
                qf2 = ((2.0 * k - 1.0) * x * qf1 - (k - 1.0) * qf0) / k
                qm[k+1] = qf2
                qf0 = qf1
                qf1 = qf2
            end
        end

        qg0 = q01
        qg1 = q11
        for k in 2:m
            qm1 = -2.0 * (k - 1.0) / xq * x * qg1 - ls * k * (3.0 - k) * qg0
            qg0 = qg1
            qg1 = qm1
        end

        if m == 0
            qm1 = q01
        elseif m == 1
            qm1 = q11
        end

        qm[2] = qm1

        if (m == 1) && (n > 1)
            qh0 = q10
            qh1 = q11
            for k in 2:n
                qh2 = ((2.0 * k - 1.0) * x * qh1 - k * qh0) / (k - 1.0)
                qm[k+1] = qh2
                qh0 = qh1
                qh1 = qh2
            end
        elseif m >= 2
            qg0 = q00
            qg1 = q01
            qh0 = q10
            qh1 = q11
            qmk = 0.0
            for l in 2:n
                q0l = ((2.0 * l - 1.0) * x * qg1 - (l - 1.0) * qg0) / l
                q1l = ((2.0 * l - 1.0) * x * qh1 - l * qh0) / (l - 1.0)
                qf0 = q0l
                qf1 = q1l
                for k in 2:m
                    qmk = -2.0 * (k - 1.0) / xq * x * qf1 - ls * (k + l - 1.0) * (l + 2.0 - k) * qf0
                    qf0 = qf1
                    qf1 = qmk
                end
                qm[l+1] = qmk
                qg0 = qg1
                qg1 = q0l
                qh0 = qh1
                qh1 = q1l
            end
        end
    else
        if abs(x) > 1.1
            km = 40 + m + n
        else
            km = trunc(Int, -1.0 - 1.8 * log(x - 1.0)) * (40 + m + n)
        end
        qf2 = 0.0
        qf1 = 1.0
        for k in km:-1:0
            qf0 = ((2.0 * k + 3.0) * x * qf1 - (k + 2.0 - m) * qf2) / (k + m + 1.0)
            if k <= n
                qm[k+1] = qf0
            end
            qf2 = qf1
            qf1 = qf0
        end
        for k in 0:n
            qm[k+1] *= qm0 / qf0
        end
    end

    if abs(x) < 1.0
        for k in 0:n
            qm[k+1] *= (-1) ^ m
        end
    end

    qd[1] = ((1.0 - m) * qm[2] - x * qm[1]) / (x * x - 1.0)
    for k in 1:n
        qd[k+1] = (k * x * qm[k+1] - (k + m) * qm[k]) / (x * x - 1.0)
    end
end


"""
Compute the associated Legendre function
Pmv(x) with an integer order and an arbitrary
nonnegative degree v

Input:
- v   --- Degree of Pmv(x),     v >= 0
- m   --- Order of Pmv(x)
- x   --- Argument of Pm(x),    -1 ≤ x ≤ 1

Output:
- PMV --- Pmv(x)

Routine called:
- PSI for computing Psi function
"""
function lpmv0(v::T, m::Int, x::T)::T where {T<:AbstractFloat}
    @assert v >= 0
    @assert abs(x) <= 1
    el = SF_CONST_EULER
    _EPS = 1e-14

    nv = trunc(Int, v)
    v0 = v - nv
    if x == -1.0 && v != nv
        if m == 0
            return -1.0e300
        elseif m != 0
            return 1.0e300
        end
    end

    c0 = 1.0
    if m != 0
        rg = v * (v + m)
        for j in 1:(m-1)
            rg *= v * v - j * j
        end
        xq = sqrt(1.0 - x * x)
        r0 = 1.0
        for j in 1:m
            r0 = 0.5 * r0 * xq / j
        end
        c0 = r0 * rg
    end

    if v0 == 0.0
        # DLMF 14.3.4, 14.7.17, 15.2.4
        pmv = 1.0
        r = 1.0
        for k in 1:(nv-m)
            r = 0.5 * r * (-nv + m + k - 1.0) * (nv + m + k) / (k * (k + m)) * (1.0 + x)
            pmv += r
        end
        return (-1) ^ nv * c0 * pmv
    end

    if x >= -0.35
        # DLMF 14.3.4, 15.2.1
        pmv = 1.0
        r = 1.0
        for k in 1:100
            r = 0.5 * r * (-v + m + k - 1.0) * (v + m + k) / (k * (m + k)) * (1.0 - x)
            pmv += r
            if k > 12 && abs(r / pmv) < _EPS
                break
            end
        end
        return (-1) ^ m * c0 * pmv
    end

    # DLMF 14.3.5, 15.8.10
    vs = sin(v * pi) / pi
    pv0 = 0.0
    if m != 0
        qr = sqrt((1.0 - x) / (1.0 + x))
        r2 = 1.0
        for j in 1:m
            r2 *= qr * j
        end
        s0 = 1.0
        r1 = 1.0
        for k in 1:(m-1)
            r1 = 0.5 * r1 * (-v + k - 1) * (v + k) / (k * (k - m)) * (1.0 + x)
            s0 += r1
        end
        pv0 = -vs * r2 / m * s0
    end

    pa = 2.0 * (psi(v) + el) + pi / tan(pi * v) + 1.0 / v
    s1 = 0.0
    for j in 1:m
        s1 += (j * j + v * v) / (j * (j * j - v * v))
    end
    pmv = pa + s1 - 1.0 / (m - v) + log(0.5 * (1.0 + x))
    r = 1.0
    for k in 1:100
        r = 0.5 * r * (-v + m + k - 1.0) * (v + m + k) / (k * (k + m)) * (1.0 + x)
        s = 0.0
        for j in 1:m
            s += ((k + j) * (k + j) + v * v) / ((k + j) * ((k + j) * (k + j) - v * v))
        end

        s2 = 0.0
        for j in 1:k
            s2 += 1.0 / (j * (j * j - v * v))
        end

        pss = pa + s + 2.0 * v * v * s2 - 1.0 / (m + k - v) + log(0.5 * (1.0 + x))
        r2 = pss * r
        pmv += r2
        if abs(r2 / pmv) < _EPS
            break
        end
    end

    return pv0 + pmv * vs * c0
end

"""
Compute the associated Legendre function
Pmv(x) with an integer order and an arbitrary
degree v, using recursion for large degrees

Input:
- v   --- Degree of Pmv(x)
- m   --- Order of Pmv(x)
- x   --- Argument of Pm(x), -1 ≤ x ≤ 1

Output:
- PMV --- Pmv(x)

Routine called:
- LPMV0
- GAMMA2
"""
function lpmv(v::T, m::Int, x::T)::T where {T<:AbstractFloat}
    @assert abs(x) <= 1
    vx = v
    mx = m
    pmv = 0.0

    if (x == -1.0) && (v != trunc(Int, v))
        if m == 0
            pmv = -Inf
        else
            pmv = Inf
        end
        return pmv
    end

    # DLMF 14.9.5
    if v < 0
        vx = -vx - 1.0
    end
    neg_m = 0
    if m < 0
        if ((vx + m + 1) > 0) || (vx != trunc(Int, vx))
            neg_m = 1
            mx = -m
        else
            # We don't handle cases where DLMF 14.9.3 doesn't help
            return NaN
        end
    end

    nv = trunc(Int, vx)
    v0 = vx - nv
    if (nv > 2) && (nv > mx)
        # Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3
        p0 = lpmv0(v0 + mx, mx, x)
        p1 = lpmv0(v0 + mx + 1, mx, x)
        pmv = p1
        for j in (mx + 2):nv
            pmv = ((2 * (v0 + j) - 1) * x * p1 - (v0 + j - 1 + mx) * p0) / (v0 + j - mx)
            p0 = p1
            p1 = pmv
        end
    else
        pmv = lpmv0(vx, mx, x)
    end

    if (neg_m != 0) && (abs(pmv) < 1.e300)
        # DLMF 14.9.3
        g1 = gamma2(vx - mx + 1)
        g2 = gamma2(vx + mx + 1)
        pmv *= g1 / g2 * (-1) ^ mx
    end

    return pmv
end
