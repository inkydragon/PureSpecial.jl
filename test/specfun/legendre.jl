# SPDX-License-Identifier: MIT

@testset "_lpmv0" begin
    test_v = Float64[
        1.1, 3.14,
        1:5...,
    ]
    test_m = Int[
        0:5...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        # `!(x >= -0.35)`
        -0.4,
    ]

    for v in test_v,
        x in test_x,
        m in test_m
        @assert abs(x) <= 1
        @testset "_lpmv0(v=$v, m=$m, x=$x)" begin
            p_ref = _lpmv0(v, m, x)
            p = Specfun.lpmv0(v, m, x)
            @test isapprox(p_ref, p)
        end   
    end
end

@testset "lpmv" begin
    test_v = Float64[
        -1:5...,
        1.1, 3.14,
        1:5...,
    ]
    test_m = Int[
        -4:4...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        # `!(x >= -0.35)`
        -0.4,
    ]

    for v in test_v,
        x in test_x,
        m in test_m
        @assert abs(x) <= 1
        @testset "lpmv(v=$v, m=$m, x=$x)" begin
            p_ref = _lpmv(v, m, x)
            p = Specfun.lpmv(v, m, x)
            @test isapprox(p_ref, p; nans=true)
        end   
    end
end
