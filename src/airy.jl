# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Airy functions
airy(z)
	Airy functions and their derivatives.
	+ airy_wrap -> cephes_airy
	+ cairy_wrap -> amos_airy,amos_biry

airye(z)
	Exponentially scaled Airy functions and their derivatives.
	+ cairy_wrap_e -> amos_airy,amos_biry
	+ cairy_wrap_e_real -> amos_airy,amos_biry

ai_zeros(nt)
	+ _specfun.airyzo(nt, 1) -> specfun_airyzo

bi_zeros(nt)
	+ _specfun.airyzo(nt, 2) -> specfun_airyzo

itairy(x)
    Integrals of Airy functions
    + itairy_wrap -> specfun_itairy
"""

export ai_zeros, bi_zeros


"""
	ai_zeros(nt)

Compute the first `nt` zeros of Airy functions `Ai(x)` and `Ai'(x)`, `a` and `a'`,
and the associated values of `Ai(a')` and `Ai'(a)`.
"""
function ai_zeros(nt)
	if nt < 0
		throw(DomainError(nt, "`nt` must be a positive integer (nt >= 0)"))
	end

	# for Ai(x) and Ai'(x)
	kf = 1
	a, b, c, d = zeros(nt), zeros(nt), zeros(nt), zeros(nt)

	Specfun.airyzo!(nt, kf, a, b, c, d)

	a, b, c, d
end

"""
	bi_zeros(nt)

Compute the first `nt` zeros of Airy functions `Bi(x)` and `Bi'(x)`, `b` and `b'`,
and the associated values of `Bi(b')` and `Bi'(b)`.
"""
function bi_zeros(nt)
	if nt < 0
		throw(DomainError(nt, "`nt` must be a positive integer (nt >= 0)"))
	end

	# for Bi(x) and Bi'(x)
	kf = 2
	a, b, c, d = zeros(nt), zeros(nt), zeros(nt), zeros(nt)

	Specfun.airyzo!(nt, kf, a, b, c, d)

	a, b, c, d
end
