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
- 🚧 lpmns:    Associated legendre function for a given order, Pmn(x), Pmn'(x) 
- lqmn:     Associated legendre function, Qmn(x), Qmn'(x)
- clqmn:    Associated legendre function, Qmn(z), Qmn'(z)
- 🚧 lqmns:    Associated legendre function for a given order, Qmn(x), Qmn'(x)
- 🚧 lpmv:     Associated legendre function, Pmv(x), |x| <= 1, v isa Real, v >= 0, 
"""

"""
Compute associated Legendre functions Pmn(x)
and Pmn'(x) for a given order

Input :
- x --- Argument of Pmn(x)
- m --- Order of Pmn(x),  m = 0,1,2,...,n
- n --- Degree of Pmn(x), n = 0,1,2,...,N

Output:
- PM(n) --- Pmn(x)
- PD(n) --- Pmn'(x)
"""
function lpmns()

end

"""
Compute associated Legendre functions Qmn(x)
and Qmn'(x) for a given order

Input :  
- x --- Argument of Qmn(x)
- m --- Order of Qmn(x),  m = 0,1,2,...
- n --- Degree of Qmn(x), n = 0,1,2,...

Output:  
- QM(n) --- Qmn(x)
- QD(n) --- Qmn'(x)
"""
function lqmns()

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

Input :  
- x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
- m   --- Order of Pmv(x)
- v   --- Degree of Pmv(x)

Output:  
- PMV --- Pmv(x)

Routine called:  
- LPMV0
- GAMMA2
"""
function lpmv()

end
