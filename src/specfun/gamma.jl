# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Gamma Functions"""
#=
- ✅ gamma2
- ✅ lgama
- ✅ cgama

- ✅ beta
- ✅ psi
- ✅ cpsi

- ✅ incog
- ✅ incob

- ✅ gam0
- ✅ gaih
=#

"""
Ref
- DLMF 5.7.1,5.7.2
- Wrench, J. W. (1968).
    Concerning two series for the gamma function.
    Mathematics of Computation, 22(103), 617-626.
    https://www.jstor.org/stable/2004538

TODO: The coefficients are not consistent with those in the paper
"""
const _GAM0_G = NTuple{25, Float64}((
#= 01 =#
    1.0e0, 0.57721_56649_01532_9e0,
    #                    9
    -0.65587_80715_20253_8e0, -0.4200_26350_34095_2e-1,
     0.16653_86113_82291_5e0, -0.4219_77345_55544_3e-1,
    #                6_973                     1
    -0.962_19715_27877_0e-2, 0.721_89432_46663_0e-2,
    #                                          95098
    -0.116_51675_91859_1e-2, -0.21_52416_74114_9e-3,
#= 11 =#
    #                 1162                   78824
     0.12_80502_82388_2e-3, -0.2_01348_54780_7e-4,
    -0.12504934821e-5, 0.11330272320e-5, 
    -0.2056338417e-6, 0.61160950e-8, 
     0.50020075e-8, -0.11812746e-8,
     0.1043427e-9, 0.77823e-11,
#= 21 =#
    -0.36968e-11, 0.51e-12, 
    -0.206e-13, -0.54e-14,
    0.14e-14
))

"""
    gam0(x)

Compute gamma function Г(x) for small `x`

Input
x  --- Argument of Г(x)  ( |x| ≤ 1 )

Output
ga --- Г(x)
"""
function gam0(x::Float64)
    @assert abs(x) <= 1

    gr = _GAM0_G[end]

    for k = 24:-1:1
        gr = gr * x + _GAM0_G[k]
    end

    return 1.0 / (gr * x)
end

const _GAMMA2_G = NTuple{26, Float64}((
    _GAM0_G...,
    0.1e-15,
))

"""
    gamma2(x::Float64)

Compute gamma function `Г(x)`.

Input
- `x`  --- Argument of Г(x)
        ( x is not equal to 0,-1,-2,… )

Output
- Г(x)
"""
function gamma2(x::Float64)
    ga = NaN
    if isinteger(x)
        if x > 0.0
            # When x is positive integer
            #   DLMF 5.4.1:  Г(n+1) = n!
            # CoSF (3.1.5)
            ga = 1.0
            m1 = trunc(Int64, x) - 1
            for k = 2:m1
                ga *= k
            end
        else
            # Set Г(x) = Inf, when x = -n
            ga = SF_INF300
        end
    else
        r = 1.0
        z = 0.0
        if abs(x) > 1.0
            # When |x| > 1, DLMF 5.5.1:
            #   Г(z+1) = z*Г(z)
            # CoSF (3.1.9)
            z = abs(x)
            m = trunc(Int64, z)
            for k = 1:m
                r *= (z - k)
            end
            z -= m
        else
            # When |x| <= 1
            z = x
        end

        # Calculate 1/Г(x), DLMF 5.7.1:
        # CoSF (3.1.15)
        gr = _GAMMA2_G[26]
        @inbounds for k = 25:-1:1
            gr = gr*z + _GAMMA2_G[k]
        end
        ga = 1.0 / (gr * z)

        if abs(x) > 1.0
            # When |x| > 1, DLMF 5.5.1:
            # CoSF (3.1.9)
            ga *= r
            if x < 0.0
                # When x < 0, DLMF 5.5.3:
                # CoSF (3.1.10)
                ga = -pi / (x * ga * sin(pi * x))
            end
        end
    end

    return ga
end

"""
Coefficients for the series expansion
"""
const _LGAMA_A = NTuple{10, Float64}((
    8.333333333333333e-02, -2.777777777777778e-03,
    7.936507936507937e-04, -5.952380952380952e-04,
    8.417508417508418e-04, -1.917526917526918e-03,
    6.410256410256410e-03, -2.955065359477124e-02,
    1.796443723688307e-01, -1.392432216905900e+00
))

"""
    lgama(kf::Int, x::Float64)

Compute gamma function `ln[Γ(x)]` or `Γ(x)`.

Input:
- `x`  --- Argument of Γ(x) ( x > 0 )
- `kf` --- Function code
    - kf=0 for `ln[Γ(x)]`
    - kf=1 for `Γ(x)`

Output:
- `ln[Γ(x)]` or `Γ(x)`
"""
function lgama(kf::Int, x::Float64)
    @assert kf in [1, 0]
    @assert x > 0

    x0 = x
    n = 0
    if x == 1.0 || x == 2.0
        # Special case
        gl = 0.0
        return ifelse(kf == 1, exp(gl), gl)
    elseif x <= 7.0
        # When x <= 7, add an integer to x,
        #   such that x+n > 7
        n = trunc(Int, 7 - x)
        x0 = x + n
    end

    # Calculate ln[Γ(x)] using CoSF (3.1.16)
    # DLMF 5.11.1: Asymptotic Expansions
    x2 = 1.0 / (x0 * x0)
    xp = 6.283185307179586477
    gl0 = _LGAMA_A[10]
    for k in 9:-1:1
        gl0 = gl0 * x2 + _LGAMA_A[k]
    end
    gl = gl0 / x0 + 0.5 * log(xp) + (x0 - 0.5) * log(x0) - x0

    if x <= 7.0
        # Calculate ln[Γ(x)] using CoSF (3.1.9)
        # DLMF 5.5.1: Recurrence Relations
        for k in 1:n
            gl -= log(x0 - 1.0)
            x0 -= 1.0
        end
    end

    return ifelse(kf == 1, exp(gl), gl)
end

"""
    gaih(x::Float64)

Compute gamma function Г(x)

Input
x  --- Argument of Г(x), x = n/2, n=1,2,3,…

Output
ga --- Г(x)
"""
function gaih(x::Float64)
    @assert x > 0

    ga = NaN
    if x > 0.0
        if isinteger(x)
            m1 = trunc(Int64, x - 1.0)
            ga = 1.0
            for k = 2:m1
                ga *= k
            end
        elseif isinteger(x + 0.5)
            m = trunc(Int64, x)
            ga = sqrt(pi)
            for k = 1:m
                ga *= 0.5 * (2.0 * k - 1.0)
            end
        end
    end

    return ga
end
gaih(x::Integer) = gaih(float(x))

"""
    cgama(z::Complex{Float64}, kf::Int)

Compute complex gamma function `ln[Г(z)]` or `Г(z)`.

Input
- `z`  --- Complex argument
- `kf` --- Function code
    - `kf=1` for `Г(z)`
    - `kf=0` for `ln[Г(z)]`

Output
- `ln[Г(z)]` or `Г(z)`
"""
function cgama(z::Complex{Float64}, kf::Int)
    @assert kf in [0, 1] "Only accept `kf=0` or `kf=1`"

    x = real(z)
    y = imag(z)

    # Check if z is a negative real integer
    if y == 0.0 && x <= 0.0 && isinteger(x)
        # Inf
        return complex(SF_INF300)
    end

    x1 = x
    y1 = y
    if x < 0.0
        # When Re(z) < 0, let z = -z
        x = -x
        y = -y
        z = -z
        @assert z === complex(x, y)
    else
        x1 = x
        y1 = 0.0
    end

    x0 = x
    na = 0
    if x <= 7.0
        # When x <= 7, add an integer to x,
        #   such that x+n > 7
        na = trunc(Int64, 7 - x)
        x0 += na
    end

    # Calculate ln[Γ(z+n)] using CoSF (3.1.16)
    # DLMF 5.11.1: Asymptotic Expansions
    az0 = sqrt(x0*x0 + y*y)
    th = atan(y / x0)
    gr = (x0 - 0.5)*log(az0) - th*y - x0 + 0.5*log(2.0 * pi)
    gi = th * (x0 - 0.5) + y*log(az0) - y
    for k in 1:10
        t = az0 ^ (1 - 2*k)
        gr += _LGAMA_A[k] * t * cos((2.0*k - 1.0) * th)
        gi -= _LGAMA_A[k] * t * sin((2.0*k - 1.0) * th)
    end

    if x <= 7.0
        # Calculate ln[Γ(z)] using CoSF (3.1.9)
        # DLMF 5.5.1: Recurrence Relations
        gr1 = 0.0
        gi1 = 0.0
        for j in 0:(na-1)
            gr1 += 0.5 * log((x + j)^2 + y^2)
            gi1 += atan(y / (x + j))
        end
        gr -= gr1
        gi -= gi1
    end

    if x1 < 0.0
        # Apply CoSF (3.1.10) when Re(z) < 0
        # DLMF 5.5.3: Reflection Relations
        az0 = sqrt(x*x + y*y)
        th1 = atan(y / x)
        sr = -sin(pi*x) * cosh(pi*y)
        si = -cos(pi*x) * sinh(pi*y)
        az1 = sqrt(sr*sr + si*si)
        th2 = atan(si / sr)
        if sr < 0.0
            th2 += pi
        end
        gr = log(pi / (az0 * az1)) - gr
        gi = -th1 - th2 - gi
        z = complex(x1, y1)
    end

    if kf == 1
        # Calculate Г(z)
        # TODO-OPT: return exp(gr) * cis(gi)
        g0 = exp(gr)
        gr = g0 * cos(gi)
        gi = g0 * sin(gi)
    end
    # else:  ln[ Г(z) ]

    return complex(gr, gi)
end
cgama(f::Float64, kf) = cgama(complex(f), kf)

"""
Compute beta function `B(p, q)`.

Input :
- `p`  --- Parameter  ( `p > 0` )
- `q`  --- Parameter  ( `q > 0` )

Output:
- `B(p, q)`

Routine called:
- [`Specfun.gamma2`](@ref) for computing `Γ(x)`
"""
function beta(p::Float64, q::Float64)
    @assert p > 0
    @assert q > 0

    gp = gamma2(p)
    gq = gamma2(q)
    ppq = p + q
    gpq = gamma2(ppq)
    # Calculate B(p,q) using CoSF (3.2.6)
    # DLMF 5.12.1:  B(a,b) = Γ(a)*Γ(b)/Γ(a+b)
    bt = gp * gq / gpq

    return bt
end

"""
Compute the incomplete gamma function
`r(a,x)`, `Γ(a,x)` and `P(a,x)`

Input :
- `a`   --- Parameter ( a ≤ 170 )
- `x`   --- Argument

Output: `(GIN, GIM, GIP, ISFER)`:
- `GIN` --- r(a,x)
- `GIM` --- Γ(a,x)
- `GIP` --- P(a,x)
- `ISFER` --- Error flag

Routine called:
- [`Specfun.gamma2`](@ref)
"""
function incog(a::Float64, x::Float64)
    gin, gim, gip, isfer = NaN, NaN, NaN, 0

    xam = -x + a * log(x)
    if xam > 700.0 || a > 170.0
        isfer = 6
        return gin, gim, gip, isfer
    end

    if x == 0.0
        gin = 0.0
        ga = gamma2(a)
        gim = ga
        gip = 0.0
    elseif x <= (1.0 + a)
        s = 1.0 / a
        r = s
        for k in 1:60
            r = r * x / (a + k)
            s = s + r
            if abs(r / s) < SF_EPS15
                break
            end
        end
        gin = exp(xam) * s
        ga = gamma2(a)
        gip = gin / ga
        gim = ga - gin
    elseif x > (1.0 + a)
        t0 = 0.0
        for k in 60:-1:1
            t0 = (k - a) / (1.0 + k / (x + t0))
        end
        gim = exp(xam) / (x + t0)
        ga = gamma2(a)
        gin = ga - gim
        gip = 1.0 - gim / ga
    end

    return gin, gim, gip, isfer
end

"""
Compute the incomplete beta function Ix(a,b)

Input :
- `a` --- Parameter
- `b` --- Parameter
- `x` --- Argument ( 0 ≤ x ≤ 1 )

Output:
- `BIX` --- Ix(a,b)

Routine called:
- [`Specfun.beta`](@ref) for computing beta function B(p,q)
"""
function incob(a::Float64, b::Float64, x::Float64)
    @assert 0 <= x <= 1
    dk = zeros(Float64, 51)
    fk = zeros(Float64, 51)

    s0 = (a + 1.0) / (a + b + 2.0)
    bt = beta(a, b)
    bix = NaN
    if x <= s0
        for k in 1:20
            dk[2*k] = k * (b - k) * x / ((a + 2.0*k - 1.0) * (a + 2.0*k))
        end
        for k in 0:20
            dk[2*k + 1] = -(a + k) * (a + b + k) * x / ((a + 2.0*k) * (a + 2.0*k + 1.0))
        end
        t1 = 0.0
        for k in 20:-1:1
            t1 = dk[k] / (1.0 + t1)
        end
        ta = 1.0 / (1.0 + t1)
        bix = x^a * (1.0 - x)^b / (a * bt) * ta
    else
        for k in 1:20
            fk[2*k] = k * (a - k) * (1.0 - x) / ((b + 2.0*k - 1.0) * (b + 2.0*k))
        end
        for k in 0:20
            fk[2*k + 1] = -(b + k) * (a + b + k) * (1.0 - x) / ((b + 2.0*k) * (b + 2.0*k + 1.0))
        end
        t2 = 0.0
        for k in 20:-1:1
            t2 = fk[k] / (1.0 + t2)
        end
        tb = 1.0 / (1.0 + t2)
        bix = 1.0 - x^a * (1.0 - x)^b / (b * bt) * tb
    end

    return bix
end


const _PSI_A = NTuple{8, Float64}((
    -0.8333333333333e-01,       0.83333333333333333e-02,
    -0.39682539682539683e-02,   0.41666666666666667e-02,
    -0.75757575757575758e-02,   0.21092796092796093e-01,
    -0.83333333333333333e-01,   0.4432598039215686
)) # _PSI_A

"""
Compute Psi function `Ψ(x)` (Digamma Function).

Input:
- `x`  --- Argument of `psi(x)`

Output:
- `psi(x)`
"""
function psi(x::T) where {T<:AbstractFloat}
    _EL = SF_EULER_GAMMA
    _2LOG2 = 1.386294361119891
    @assert isapprox(2*log(2), _2LOG2)

    if (x == trunc(Int, x)) && (x <= 0.0)
        return T(SF_INF300)
    end

    xa = abs(x)
    s = 0.0
    if xa == trunc(Int, xa)
        # Use CoSF (3.3.7) for |x| = n
        n = trunc(Int, xa)
        for k in 1:(n-1)
            s += 1.0 / k
        end
        ps = -_EL + s
    elseif (xa + 0.5) == trunc(Int, xa + 0.5)
        # Use CoSF (3.3.6) for |x| = n+1/2
        n = trunc(Int, xa - 0.5)
        for k in 1:n
            s += 1.0 / (2.0 * k - 1.0)
        end
        ps = -_EL + 2.0 * s - _2LOG2
    else
        if xa < 10.0
            # When |x| < 10, add an integer to |x|,
            #   such that |x| > 10
            n = 10 - trunc(Int, xa)
            for k in 0:(n-1)
                # Caculate the summation in CoSF (3.3.8)
                s += 1.0 / (xa + k)
            end
            xa += n
        end

        # Calculate psi(|x|+n) using CoSF (3.3.14)
        #   Asymptotic Expansions, with 8 terms
        a = _PSI_A
        x2 = 1.0 / (xa * xa)
        ps = log(xa) - 0.5 / xa
        ps += x2 *
            (((((((a[8]*x2+a[7])*x2+a[6])*x2+a[5])*x2+a[4])*x2+a[3])*x2+a[2])*x2+a[1])
        # Calaculate psi(|x|) using CoSF (3.3.8)
        #   psi(x) = psi(|x|+n) - sum
        ps -= s
    end

    if x < 0.0
        ps -= pi * cos(pi * x) / sin(pi * x) + 1.0 / x
    end

    return T(ps)
end

"""
    cpsi(z::ComplexF64)

Compute complex psi function `Ψ(z)` (Digamma Function).

Input :
- `z`

Output:
- `psi(z)`
"""
function cpsi(z::ComplexF64)
    psr, psi = NaN, NaN
    x, y = real(z), imag(z)

    if y == 0.0 && x == trunc(Int, x) && x <= 0.0
        return complex(1.0e300, 0.0)
    end

    x1 = x
    y1 = y
    if x < 0.0
        # When Re(z) < 0, let z = -z
        x = -x
        y = -y
    end
    x0 = x
    n = 0
    if x < 8.0
        # When x < 8, add an integer to x,
        #   such that x+n > 8
        n = 8 - trunc(Int, x)
        x0 += n
    end

    th = 0.0
    if x0 == 0.0 && y != 0.0
        # TODO: Unreachable
        th = 0.5 * pi
    elseif x0 != 0.0
        th = atan(y / x0)
    end

    # Calculate psi(z+n) using CoSF (3.3.14)
    z2 = x0*x0 + y*y
    z0 = sqrt(z2)
    psr = log(z0) - 0.5 * x0 / z2
    psi = th + 0.5 * y / z2
    for k in 1:8
        psr += _PSI_A[k] * z2^(-k) * cos(2.0 * k * th)
        psi -= _PSI_A[k] * z2^(-k) * sin(2.0 * k * th)
    end

    if x < 8.0
        # Calculate psi(z) using CoSF (3.3.14)
        rr = 0.0
        ri = 0.0
        for k in 1:n
            rr += (x0 - k) / ((x0 - k)^2 + y*y)
            ri += y / ((x0 - k)^2 + y*y)
        end
        psr -= rr
        psi += ri
    end

    if x1 < 0.0
        # Apply CoSF (3.3.11) whne Re(z) < 0
        tn = tan(pi * x)
        tm = tanh(pi * y)
        ct2 = tn*tn + tm*tm
        psr = psr + x / (x*x + y*y) + pi * (tn - tn * tm*tm) / ct2
        psi = psi - y / (x*x + y*y) - pi * tm * (1.0 + tn*tn) / ct2
        x = x1
        y = y1
    end

    return complex(psr, psi)
end
