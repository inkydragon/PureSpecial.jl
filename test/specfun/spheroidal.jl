# SPDX-License-Identifier: MIT

"""Spheroidal Wave Functions
- specfun_segv
- specfun_rswfp
- specfun_aswfa
- specfun_rswfo
"""

@testset "segv" begin
    # TODO: test branch: `if k != 1 && h[k] < h[k-1]`
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 20),
        (50, 100),
        (100, 200),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9, eps(), 1e-10,
        1:10...,
        rand(10)...,
    ]

    for (m, n) in test_mn,
        c in test_c,
        kd in [1, -1]
        #
        max_cv_len = n - m + 1
        @assert 0 <= max_cv_len <= 200
        icm = (n - m + 2) ÷ 2
        max_eg_nz_len = icm * 2
        eg_ref = zeros(Float64, 200)
        eg_res = zeros(Float64, 200)
        #
        cv_ref, _ = _segv!(m, n, c, kd, eg_ref)
        cv_res, _ = Specfun.segv(m, n, c, kd, eg_res)
        @testset "_segv(m=$m, n=$n, c=$c, kd=$kd)" begin
            # Result
            @test isapprox(cv_ref, cv_res)
            @test isapprox(eg_ref[1:max_cv_len], eg_res[1:max_cv_len])
            # temp value
            @test isapprox(eg_ref[max_cv_len+1:max_eg_nz_len], eg_res[max_cv_len+1:max_eg_nz_len])
            # zeros
            @test iszero(eg_ref[max_eg_nz_len+1:end])
            @test iszero(eg_res[max_eg_nz_len+1:end])
        end
    end
end
