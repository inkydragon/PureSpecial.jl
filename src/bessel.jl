# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Bessel functions
jv(v, z)
    Bessel function of the first kind of real order and complex argument.
    + cbesj_wrap_real -> cbesj_wrap,cephes_jv
    + cbesj_wrap -> amos_besj,amos_besy
        + cbesj_wrap_e

jve(v, z)
    Exponentially scaled Bessel function of the first kind of order v.
    + cbesj_wrap_e_real -> cbesj_wrap_e
    + cbesj_wrap_e -> amos_besj,amos_besy

yn(n, x)
    Bessel function of the second kind of integer order and real argument.
    + cephes_yn

yv(v, z)
    Bessel function of the second kind of real order and complex argument.
    + cbesy_wrap_real -> cbesy_wrap,cephes_yv
    + cbesy_wrap -> amos_besy,amos_besj

yve(v, z)
    Exponentially scaled Bessel function of the second kind of real order.
    + cbesy_wrap_e_real -> cbesy_wrap_e
    + cbesy_wrap_e -> amos_besy,amos_besj

kn(n, x)
    Modified Bessel function of the second kind of integer order n
    + cbesk_wrap_real_int -> cbesk_wrap_real -> cbesk_wrap

kv(v, z)
    Modified Bessel function of the second kind of real order v
    + cbesk_wrap_real
    + cbesk_wrap

kve(v, z)
    Exponentially scaled modified Bessel function of the second kind.
    + cbesk_wrap_e_real
    + cbesk_wrap_e

iv(v, z)
    Modified Bessel function of the first kind of real order.
    + cbesi_wrap
    + cephes_iv

ive(v, z)
    Exponentially scaled modified Bessel function of the first kind.
    + cbesi_wrap_e_real
    + cbesi_wrap_e
        
hankel1(v, z)
    Hankel function of the first kind
    + cbesh_wrap1

hankel1e(v, z)
    Exponentially scaled Hankel function of the first kind
    + cbesh_wrap1_e

hankel2(v, z)
    Hankel function of the second kind
    + cbesh_wrap2

hankel2e(v, z)
    Exponentially scaled Hankel function of the second kind
    + cbesh_wrap2_e

wright_bessel(a, b, x)
    Wright's generalized Bessel function.
    + wright_bessel_scalar


lmbda(v, x)
Jahnke-Emden Lambda function, Lambdav(x).

"""

