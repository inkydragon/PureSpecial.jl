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
Euler's Constant γ
"""
const SF_CONST_EULER = 0.5772156649015329

const _PSI_A = NTuple{8, Float64}((
    -0.8333333333333e-01,
    0.83333333333333333e-02,
    -0.39682539682539683e-02,
    0.41666666666666667e-02,
    -0.75757575757575758e-02,
    0.21092796092796093e-01,
    -0.83333333333333333e-01,
    0.4432598039215686
)) # _PSI_A

"""
Compute Psi function

Input :  
- x  --- Argument of psi(x)

Output:  
- PS --- psi(x)

TODO: move to gamma.jl
"""
function psi(x::T) where {T<:AbstractFloat}
    @assert isapprox(pi, 3.141592653589793)
    _EL = SF_CONST_EULER
    _2LOG2 = 1.386294361119891
    @assert isapprox(2*log(2), _2LOG2)

    if (x == trunc(Int, x)) && (x <= 0.0)
        return T(1e300)
    end

    xa = abs(x)
    s = 0.0
    if xa == trunc(Int, xa)
        n = trunc(Int, xa)
        for k in 1:(n-1)
            s += 1.0 / k
        end
        ps = -_EL + s
    elseif (xa + 0.5) == trunc(Int, xa + 0.5)
        n = trunc(Int, xa - 0.5)
        for k in 1:n
            s += 1.0 / (2.0 * k - 1.0)
        end
        ps = -_EL + 2.0 * s - _2LOG2
    else
        if xa < 10.0
            n = 10 - trunc(Int, xa)
            for k in 0:(n-1)
                s += 1.0 / (xa + k)
            end
            xa += n
        end
        a = _PSI_A
        x2 = 1.0 / (xa * xa)
        ps = log(xa) - 0.5 / xa
        ps += x2 *
            (((((((a[8]*x2+a[7])*x2+a[6])*x2+a[5])*x2+a[4])*x2+a[3])*x2+a[2])*x2+a[1])
        ps -= s
    end

    if x < 0.0
        ps -= pi * cos(pi * x) / sin(pi * x) + 1.0 / x
    end

    return T(ps)
end


"""
Compute the associated Legendre function
Pmv(x) with an integer order and an arbitrary
nonnegative degree v

Input :  
- x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
- m   --- Order of Pmv(x)
- v   --- Degree of Pmv(x)

Output:  
- PMV --- Pmv(x)

Routine called:  
- PSI for computing Psi function
"""
function lpmv0()

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
