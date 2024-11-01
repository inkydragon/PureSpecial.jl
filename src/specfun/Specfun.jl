# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Special functions in pure julia.

+ ✅ specfun_airyzo
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
+ ✅ specfun_jdzo
+ ✅ specfun_jyzo
+ ✅ specfun_klvnzo
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

## Airy functions
- `ai_zeros(nt)`: specfun_airyzo
- `bi_zeros(nt)`: specfun_airyzo
- `itairy(x)`: specfun_itairy
"""

module Specfun

# SUBROUTINE
include("gamma.jl")

include("airy.jl")
include("bessel_zeros.jl")
include("error.jl")
include("parabolic.jl")
include("kelvin.jl")
include("hyper.jl")
include("exp.jl")

end # PureSpecial