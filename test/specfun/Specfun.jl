# SPDX-License-Identifier: MIT
using PureSpecial.Specfun


"""Special functions.

+ specfun_airyzo
+ specfun_bernob
+ specfun_cerzo
+ specfun_clpmn
+ specfun_clpn
+ specfun_clqmn
+ specfun_clqn
+ specfun_cpbdn
+ specfun_cyzo
+ specfun_eulerb
+ specfun_fcoef
+ specfun_fcszo
+ specfun_jdzo
+ specfun_jyzo
+ specfun_klvnzo
+ specfun_lamn
+ specfun_lamv
+ specfun_lpmn
+ specfun_lpn
+ specfun_lqmn
+ specfun_lqnb
+ specfun_pbdv
+ specfun_rctj
+ specfun_rcty
+ specfun_sdmn
+ specfun_segv
"""

# for Test
include("warp.jl")

# SUBROUTINE
include("gamma.jl")

include("airy.jl")
include("bessel_zeros.jl")
include("error.jl")
include("parabolic.jl")
include("kelvin.jl")
include("hyper.jl")
include("exp.jl")
include("legendre.jl")
include("spheroidal.jl")
