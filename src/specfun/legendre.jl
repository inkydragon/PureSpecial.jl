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
- lpmns:    Associated legendre function for a given order, Pmn(x), Pmn'(x) 
- lqmn:     Associated legendre function, Qmn(x), Qmn'(x)
- clqmn:    Associated legendre function, Qmn(z), Qmn'(z)
- lqmns:    Associated legendre function for a given order, Qmn(x), Qmn'(x)
- lpmv:     Associated legendre function, Pmv(x), |x| <= 1, v isa Real, v >= 0, 
"""

