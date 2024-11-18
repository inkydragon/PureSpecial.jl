# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Special functions in pure julia.

Translated from `specfun.f`

```fortran
C
C       COMPUTATION OF SPECIAL FUNCTIONS
C
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
```



"""
module Specfun

include("const.jl")

include("gamma.jl")
include("exp.jl")
include("error.jl")
include("airy.jl")

include("bessel_zeros.jl")
include("kelvin.jl")
include("struve.jl")
include("parabolic.jl")
include("hyper.jl")
include("legendre.jl")

include("mathieu.jl")
include("spheroidal.jl")
include("trig_int.jl")

end # PureSpecial