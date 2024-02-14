# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Scipy: scipy.special: Zeros of Bessel functions

## jdzo
jnjnp_zeros(nt)
    Compute zeros of integer-order Bessel functions Jn and Jn'.
    _specfun.jdzo -> specfun_jdzo -> bjndd

## jyzo
jnyn_zeros(n, nt)
    Compute nt zeros of Bessel functions Jn(x), Jn'(x), Yn(x), and Yn'(x).
    _specfun.jyzo -> specfun_jyzo -> jyndd -> jynbh -> msta1,msta2

    jn_zeros(n, nt) -> jnyn_zeros(n, nt)[0]
        Compute zeros of integer-order Bessel functions Jn.

    jnp_zeros(n, nt) -> jnyn_zeros(n, nt)[1]
        Compute zeros of integer-order Bessel function derivatives Jn'.

    yn_zeros(n, nt) -> jnyn_zeros(n, nt)[2]
        Compute zeros of integer-order Bessel function Yn(x).

    ynp_zeros(n, nt) -> jnyn_zeros(n, nt)[3]
        Compute zeros of integer-order Bessel function derivatives Yn'(x).

## cyzo
y0_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y0(z), and derivative at each zero.
    _specfun.cyzo(nt, 0, kc) -> specfun_cyzo -> cy01

y1_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y1(z), and derivative at each zero.
    _specfun.cyzo(nt, 1, kc) -> ...

y1p_zeros(nt[, complex])
    Compute nt zeros of Bessel derivative Y1'(z), and value at each zero.
    _specfun.cyzo(nt, 2, kc) -> ...
"""


"""
    bjndd(x::Float64, n::Int)

Compute Bessel functions Jn(x) and their
first and second derivatives ( 0 <= n <= 100 )

## Input
- `x` ---  Argument of Jn(x)  ( x ≥ 0 )
- `n` ---  Order of Jn(x)

## Output
- `BJ(n+1)` ---  Jn(x)
- `DJ(n+1)` ---  Jn'(x)
- `FJ(n+1)` ---  Jn"(x)
"""
function bjndd(x::Float64, n::Int)
    # TODO: modify params: `bjndd!`
    @assert n in 0:100
    array_szie = 101
    bj = zeros(array_szie)
    dj = zeros(array_szie)
    fj = zeros(array_szie)

    m = 1
    for nt in 1:900
        # TODO: use a table?
        mt = 0.5 * log10(6.28 * nt) - nt * log10(1.36 * abs(x) / nt)
        if mt > 20
            m = nt
            break
        end
    end #= 10 =#

    #= 15 =#
    bs = 0.0
    f = 0.0
    f0 = 0.0
    f1 = 1e-35
    for k in m:-1:0
        f = 2.0 * (k + 1.0) * f1 / x - f0
        if k <= n
            bj[k+1] = f
        end
        if k % 2 == 0
            bs += 2.0 * f
        end
        f0 = f1
        f1 = f
    end #= 20 =#

    # TODO: use @.
    for k in 1:(n + 1)
        bj[k] /= bs - f
    end #= 25 =#

    dj[1] = -bj[2]
    fj[1] = -bj[1] - dj[1] / x

    for k in 1:n
        dj[k+1] = bj[k] - k * bj[k+1] / x
        fj[k+1] = (k*k / (x*x) - 1.0) * bj[k+1] - dj[k+1] / x
    end #= 30 =#

    return bj, dj, fj
end

"""
    jdzo!(nt::Int, zo::Vector{Float64}, n::Vector{Int}, m::Vector{Int}, p::Vector{Int})

Compute the zeros of Bessel functions Jn(x) and
Jn'(x), and arrange them in the order of their
magnitudes.

## Input
NT    --- Number of total zeros ( NT ≤ 1200 )

## Output
ZO(L) --- Value of the L-th zero of Jn(x)
          and Jn'(x)
N(L)  --- n, order of Jn(x) or Jn'(x) associated
          with the L-th zero
M(L)  --- m, serial number of the zeros of Jn(x)
          or Jn'(x) associated with the L-th zero
          ( L is the serial number of all the
            zeros of Jn(x) and Jn'(x) )
P(L)  --- 0 (TM) or 1 (TE), a code for designating the
          zeros of Jn(x)  or Jn'(x).
          In the waveguide applications, the zeros
          of Jn(x) correspond to TM modes and
          those of Jn'(x) correspond to TE modes

## Routine called
BJNDD for computing Jn(x), Jn'(x) and Jn''(x)
"""
function jdzo!(nt::Int, zo::Vector{Float64}, n::Vector{Int}, m::Vector{Int}, p::Vector{Int})
    #= Array Size const =#
    # Note: ZO and ZOC arrays are 0-indexed in specfun.f
    ZO_size = 1400 + 1
    N_size = 1400
    M_size = 1400
    P_size = 1400
    @assert nt in 1:1200
    @assert length(zo) <= ZO_size
    @assert length(n) <= N_size
    @assert length(m) <= M_size
    @assert length(p) <= P_size

    ZOC_szie = 70 + 1
    N1_szie = 70
    M1_szie = 70
    P1_size = 70
    
    BJ_size = 101
    DJ_size = 101
    FJ_size = 101

    #= Initialize arrays =#
    p1 = zeros(Int64, P1_size)
    n1 = zeros(Float64, N1_szie)
    m1 = zeros(Float64, M1_szie)
    zoc = zeros(Float64, ZOC_szie)

    xm = 0.0
    nm, mm = 0, 0
    x = 0.0
    if nt < 600
        xm = -1.0 + 2.248485*sqrt(nt) - 0.0159382*nt + 3.208775e-4*nt^1.5
        nm = trunc(Int64, 14.5 + 0.05875*nt)
        mm = trunc(Int64, 0.02*nt) + 6
    else
        xm = 5.0 + 1.445389*sqrt(nt) + 0.01889876*nt - 2.147763e-4*nt^1.5
        nm = trunc(Int64, 27.8 + 0.0327*nt)
        mm = trunc(Int64, 0.01088*nt) + 10
    end

    L0 = 0
    for i in 1:nm
        x1 = 0.407658 + 0.4795504*sqrt(i-1.0) + 0.983618*(i-1)
        x2 = 1.99535  + 0.8333883*sqrt(i-1.0) + 0.984584*(i-1)
        L1 = 0
        for j in 1:mm
            if (i != 1) || (j != 1)
                x = x1
                while true
                    #= 10 =#
                    bj, dj, fj = bjndd(x, i)
                    x0 = x
                    x -= dj[i] / fj[i]
                    if x1 > xm
                        @goto _line15
                    end
                    if abs(x-x0) <= 1e-10
                        break
                    end
                end
            end

            #= 15 =#
            L1 += 1
            n1[L1] = i - 1
            m1[L1] = j
            if i == 1
                m1[L1] = j - 1
            end
            p1[L1] = 1
            zoc[L1+1] = x
            if i <= 15
                x1 = x + 3.057 + 0.0122*(i-1) + (1.555 + 0.41575*(i-1))/(j+1)^2.0
            else
                x1 = x + 2.918 + 0.01924*(i-1) + (6.26 + 0.13205*(i-1))/(j+1)^2.0
            end

            #= 20 =#
            @label _line15
            x = x2
            while true
                #= 25 =#
                bj, dj, fj = bjndd(x, i)
                x0 = x
                x -= bj[i] / dj[i]
                if x > xm
                    # Need to "continue;" twice hence goto is simpler
                    @goto _line30
                end
                if abs(x-x0) <= 1e-10
                    break
                end
            end
            L1 += 1
            n1[L1] = i - 1
            m1[L1] = j
            p1[L1] = 0
            zoc[L1+1] = x
            if i <= 15
                x2 = x + 3.11 + 0.0138*(i-1) + (0.04832 + 0.2804*(i-1))/(j+1)^2
            else
                x2 = x + 3.001 + 0.0105*(i-1) + (11.52 + 0.48525*(i-1))/(j+3)^2
            end
        end

        #= 30 =#
        @label _line30
        L = L0 + L1
        L2 = L
        while true
            #= 35 =#
            if L0 == 0
                for k in 1:L
                    p[k] = p1[k]
                    m[k] = m1[k]
                    n[k] = n1[k]
                    # NOTE: zo[] is 0-indexed
                    zo[k+1] = zoc[k+1]
                end #= 40 =#
                L1 = 0
            elseif L0 != 0
                if zo[L0+1] >= zoc[L1+1]
                    p[L0+L1] = p[L0]
                    m[L0+L1] = m[L0]
                    n[L0+L1] = n[L0]
                    zo[L0+L1+1] = zo[L0+1]
                    L0 -= 1
                else
                    p[L0+L1] = p1[L1]
                    m[L0+L1] = m1[L1]
                    n[L0+L1] = n1[L1]
                    zo[L0+L1+1] = zoc[L1+1]
                    L1 -= 1
                end
            end
            if L1 == 0
                break
            end
        end
        L0 = L2
    end #= 45 =#
end

"""
    jdzo(nt::Int64)
"""
function jdzo(nt::Int64)
    zo = zeros(Float64, 1401)
    n = zeros(Int64, 1400)
    m = zeros(Int64, 1400)
    p = zeros(Int64, 1400)

    jdzo!(nt, zo, n, m, p)
    
    zo[1:nt], n[1:nt], m[1:nt], p[1:nt]
end

"""
Helper function used in `msta1`, `msta2`
"""
function envj(n, x)
    0.5 * log10(6.28 * n) - n * log10(1.36 * x / n)
end

"""
    msta1(x::Float64, mp::Int)

Determine the starting point for backward
recurrence such that the magnitude of
Jn(x) at that point is about 10^(-MP)

Input
x     --- Argument of Jn(x)
MP    --- Value of magnitude

Output
MSTA1 --- Starting point
"""
function msta1(x::Float64, mp::Int)
    a0 = abs(x)
    n0 = trunc(Int64, 1.1 * a0) + 1
    f0 = envj(n0, a0) - mp

    n1 = n0 + 5
    f1 = envj(n1, a0) - mp

    nn = 0
    for _ = 1:20
        # NOTE: Sometimes `_nn` may be NaN.
        _nn = n1 - (n1 - n0) / (1.0 - f0 / f1)
        nn = trunc(Int64, _nn)
        f = envj(nn, a0) - mp
        if abs(nn - n1) < 1
            break
        end

        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
    end

    nn
end

"""
    msta2(x::Float64, n::Int64, mp::Int64)

Determine the starting point for backward
recurrence such that all Jn(x) has MP
significant digits

Input
x  --- Argument of Jn(x)
n  --- Order of Jn(x)
MP --- Significant digit

Output
MSTA2 --- Starting point
"""
function msta2(x::Float64, n::Int64, mp::Int64)
    a0 = abs(x)
    hmp = 0.5 * mp
    ejn = envj(n, a0)

    obj = 0.0
    n0 = 0
    if ejn <= hmp
        obj = mp
        n0 = trunc(Int64, 1.1 * a0) + 1
    else
        obj = hmp + ejn
        n0 = n
    end

    f0 = envj(n0, a0) - obj
    n1 = n0 + 5
    f1 = envj(n1, a0) - obj
    nn = 0
    for _ = 1:20
        _nn = n1 - (n1 - n0) / (1.0 - f0 / f1)
        nn = trunc(Int64, _nn)
        f = envj(nn, a0) - obj
        if abs(nn - n1) < 1
            break
        end

        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
    end

    nn + 10
end

"""
    jynbh!(
        x::Float64, n::Int, nmin::Int,
        bj::Vector{Float64}, by::Vector{Float64}
    )

Compute Bessel functions Jn(x), Yn(x)

Input
x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
nmin -- Lowest order computed  ( nmin ≥ 0 )

Output
BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0

Return:
NM --- Highest order computed

Routines called
MSTA1 and MSTA2 to calculate the starting
point for backward recurrence
"""
function jynbh!(
    x::Float64, n::Int, nmin::Int,
    bj::Vector{Float64}, by::Vector{Float64})
    out_len = length(nmin:n)
    @assert length(bj) >= out_len
    @assert length(by) >= out_len

    r2p = 0.63661977236758
    # Hankel expansion
    a = [-0.0703125, 0.112152099609375, -0.5725014209747314, 0.6074042001273483e+01]
    b = [0.0732421875, -0.2271080017089844, 0.1727727502584457e+01, -0.2438052969955606e+02]
    a1 = [0.1171875, -0.144195556640625, 0.6765925884246826, -0.6883914268109947e+01]
    b1 = [-0.1025390625, 0.2775764465332031, -0.1993531733751297e+01, 0.2724882731126854e+02]

    nm = n
    if x < 1.0e-100
        # TODO: use fill!
        for k = nmin:n
            # NOTE: in Fortran `bj,by` index start from 0
            bj[k - nmin + 1] = 0.0
            by[k - nmin + 1] = -1.0e+300
        end
        if nmin == 0
            bj[1] = 1.0
        end

        return nm
    end

    by0 = 0.0
    by1 = 0.0
    ky = 0
    if (x <= 300.0) || (n > trunc(Int64, 0.9 * x))
        #
        # Backward recurrence for Jn
        #
        if n == 0
            nm = 1
        end
        m = msta1(x, 200)
        if m < nm
            nm = m
        else
            m = msta2(x, nm, 15)
        end
    
        bs = 0.0
        su = 0.0
        sv = 0.0
        f2 = 0.0
        f1 = 1.0e-100
        f = 0.0
        for k = m:-1:0
            f = 2.0 * (k + 1.0) / x * f1 - f2
            if (k <= nm) && (k >= nmin)
                bj[k - nmin + 1] = f
            end
            if k == 2 * div(k, 2) && k != 0
                bs += 2.0 * f
                su += (-1)^(div(k, 2)) * f / k
            elseif k > 1
                sv += (-1)^(div(k, 2)) * k / (k * k - 1.0) * f
            end
            f2 = f1
            f1 = f
        end
        s0 = bs + f

        for k = nmin:nm
            bj[k - nmin + 1] /= s0
        end

        #
        # Estimates for Yn at start of recurrence
        #
        bj0 = f1 / s0
        bj1 = f2 / s0
        ec = log(x / 2.0) + 0.5772156649015329
        by0 = r2p * (ec * bj0 - 4.0 * su / s0)
        by1 = r2p * ((ec - 1.0) * bj1 - bj0 / x - 4.0 * sv / s0)

        if 0 >= nmin
            by[0 - nmin + 1] = by0
        end
        if 1 >= nmin
            by[1 - nmin + 1] = by1
        end
        ky = 2
    else
        #
        # Hankel expansion
        #
        t1 = x - 0.25 * pi
        p0 = 1.0
        q0 = -0.125 / x

        for k = 1:4
            p0 += a[k] * x^(-2*k)
            q0 += b[k] * x^(-2*k - 1)
        end

        cu = sqrt(r2p / x)
        bj0 = cu * (p0 * cos(t1) - q0 * sin(t1))
        by0 = cu * (p0 * sin(t1) + q0 * cos(t1))

        if 0 >= nmin
            bj[0 - nmin + 1] = bj0
            by[0 - nmin + 1] = by0
        end

        t2 = x - 0.75 * pi
        p1 = 1.0
        q1 = 0.375 / x

        for k = 1:4
            p1 += a1[k] * x^(-2*k)
            q1 += b1[k] * x^(-2*k - 1)
        end

        bj1 = cu * (p1 * cos(t2) - q1 * sin(t2))
        by1 = cu * (p1 * sin(t2) + q1 * cos(t2))

        if 1 >= nmin
            bj[1 - nmin + 1] = bj1
            by[1 - nmin + 1] = by1
        end

        for k = 2:nm
            bjk = 2.0 * (k - 1.0) / x * bj1 - bj0
            if k >= nmin
                bj[k - nmin + 1] = bjk
            end
            bj0 = bj1
            bj1 = bjk
        end
        ky = 2
    end

    #
    # Forward recurrence for Yn
    #
    for k = ky:nm
        byk = 2.0 * (k - 1.0) * by1 / x - by0

        if k >= nmin
            by[k - nmin + 1] = byk
        end
        by0 = by1
        by1 = byk
    end
    
    nm
end

"""
jyndd(n::Int, x::Float64)

compute bessel functions jn(x) and yn(x), and
their first and second derivatives

input
x   ---  argument of jn(x) and yn(x) ( x > 0 )
n   ---  order of jn(x) and yn(x)

output
bjn ---  jn(x)
djn ---  jn'(x)
fjn ---  jn"(x)
byn ---  yn(x)
dyn ---  yn'(x)
fyn ---  yn"(x)

routines called
jynbh to compute jn and yn
"""
function jyndd(x::Float64, n::Int)
    @assert x > 0.0
    bj, by = zeros(Float64, 2), zeros(Float64, 2)

    jynbh!(x, n + 1, n, bj, by)

    # Compute derivatives by differentiation formulas
    bjn = bj[1]
    byn = by[1]
    djn = -bj[2] + n * bj[1] / x
    dyn = -by[2] + n * by[1] / x
    fjn = (n * n / (x * x) - 1.0) * bjn - djn / x
    fyn = (n * n / (x * x) - 1.0) * byn - dyn / x
    
    return bjn, djn, fjn, byn, dyn, fyn
end
