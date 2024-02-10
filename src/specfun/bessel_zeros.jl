# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Scipy: scipy.special: Zeros of Bessel functions

## jdzo
jnjnp_zeros(nt)
    Compute zeros of integer-order Bessel functions Jn and Jn'.

## jyzo
jnyn_zeros(n, nt)
    Compute nt zeros of Bessel functions Jn(x), Jn'(x), Yn(x), and Yn'(x).

    jn_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel functions Jn.

    jnp_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function derivatives Jn'.

    yn_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function Yn(x).

    ynp_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function derivatives Yn'(x).

## cyzo
y0_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y0(z), and derivative at each zero.

y1_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y1(z), and derivative at each zero.

y1p_zeros(nt[, complex])
    Compute nt zeros of Bessel derivative Y1'(z), and value at each zero.
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

    # Compared to specfun.f we use a single array instead of separate
    # three arrays and use pointer arithmetic to access. Their usage
    # is pretty much one-shot hence does not complicate the code.

    # m1, n1, zoc -> 70 + 70 + 71
    mnzoc = zeros(Float64, M1_szie + N1_szie + ZOC_szie)

    # bj, dj, fj -> 101 + 101 + 101
    bdfj = zeros(Float64, BJ_size + DJ_size + FJ_size)

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
                    bdfj[1:101], bdfj[102:202], bdfj[203:303] = bjndd(x, i)
                    x0 = x
                    x -= bdfj[101+i] / bdfj[202+i]
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
            #= N1(L1) =#
            mnzoc[70 + L1] = i - 1
            #= M1(L1) =#
            mnzoc[L1] = j
            if i == 1
                mnzoc[L1] = j - 1
            end
            p1[L1] = 1
            #= ZOC[L1] =#
            mnzoc[141 + L1] = x
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
                bdfj[1:101], bdfj[102:202], bdfj[203:303] = bjndd(x, i)
                x0 = x
                #= X -= DJ(I) / FJ(I) =#
                x -= bdfj[i] / bdfj[101+i]
                if x > xm
                    # Need to "continue;" twice hence goto is simpler
                    @goto _line30
                end
                if abs(x-x0) <= 1e-10
                    break
                end
            end
            L1 += 1
            #= N1(L1) =#
            mnzoc[70 + L1] = i - 1
            #= M1(L1) =#
            mnzoc[L1] = j
            p1[L1] = 0
            #= ZOC[L1] =#
            mnzoc[141 + L1] = x
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
                    m[k] = trunc(Int64, mnzoc[k])
                    n[k] = trunc(Int64, mnzoc[70+k])
                    # NOTE: zo[] is 0-indexed
                    zo[k+1] = mnzoc[141+k]
                end #= 40 =#
                L1 = 0
            elseif L0 != 0
                if zo[L0+1] >= mnzoc[141+L1]
                    p[L0+L1] = p[L0]
                    m[L0+L1] = m[L0]
                    n[L0+L1] = n[L0]
                    zo[L0+L1+1] = zo[L0+1]
                    L0 -= 1
                else
                    p[L0+L1] = p1[L1]
                    m[L0+L1] = trunc(Int64, mnzoc[L1])
                    n[L0+L1] = trunc(Int64, mnzoc[70+L1])
                    zo[L0+L1+1] = mnzoc[141+L1]
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
