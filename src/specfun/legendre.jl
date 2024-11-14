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
