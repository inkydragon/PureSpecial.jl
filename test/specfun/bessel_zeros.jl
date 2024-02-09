# SPDX-License-Identifier: MIT

@testset "bjndd" begin
    test_x = Float64[
        0.0, -0.0,
        1:10...,
        rand(10:100, 5)...,
    ]

    for x in test_x,
        n in 1:5
        @testset "_bjndd(x=$x, n=$n)" begin
            r_bj, r_dj, r_fj = _bjndd(x, n)
            bj, dj, fj = Specfun.bjndd(x, n)

            @test isapprox(r_bj, bj)
            @test isapprox(r_dj, dj)
            @test isapprox(r_fj, fj)
        end
    end
end
