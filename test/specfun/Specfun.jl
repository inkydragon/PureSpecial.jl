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
include("spheroidal.jl")

@testset "sdmn" begin
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 11),
        (50, 100),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        # 1e-9, eps(), 1e-10,
        1:10...,
        # rand(10)..., # TODO:err
    ]
    test_cv = Float64[
        4.5,
        rand(10)...,
        1:10...,
        # test branch: `kb > 2 && if abs(f) > T(1e100)`
        930:950...,
    ]
    test_x = Float64[
        3.14,
        2:6...,
    ]

    for (m, n) in test_mn,
        c in test_c,
        cv in test_cv,
        x in test_x
        kd = 1
        #
        max_df_len = 1 + 25 + trunc(Int, 0.5*(n-m) + c)
        @assert 0 <= max_df_len <= 200
        @testset "_rmn2l(m=$m,n=$n, c=$c,x=$x, kd=$kd)" begin
            df_ref = zeros(Float64, 200)
            df = zeros(Float64, 200)
            #
            _sdmn!(m, n, c, cv, kd, df_ref)
            r2f_ref, r2d_ref, id_ref = _rmn2l(m, n, c, x, kd, df_ref)
            Specfun.sdmn!(m, n, c, cv, kd, df)
            r2f, r2d, id = Specfun.rmn2l(m, n, c, x, kd, df)
            # Result
            @test isapprox(r2f_ref, r2f; nans=true)
            @test isapprox(r2d_ref, r2d; nans=true)
            @test id_ref == id
        end
    end
end
