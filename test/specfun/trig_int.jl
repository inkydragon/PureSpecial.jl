# SPDX-License-Identifier: MIT

@testset "cisia" begin
    test_x = Float64[
        rand(10)...,
        0:64...,
    ]
    for x in test_x
        @testset "cisia($x)" begin
            ci_ref, si_ref = _cisia(x)
            ci, si = Specfun.cisia(x)
            @test isapprox(ci_ref, ci)
            @test isapprox(si_ref, si)
        end
    end
end

@testset "cisib" begin
    test_x = Float64[
        0:10...,
        rand(10)...,
    ]
    for x in test_x
        @testset "cisib($x)" begin
            ci_ref, si_ref = _cisib(x)
            ci, si = Specfun.cisib(x)
            @test isapprox(ci_ref, ci)
            @test isapprox(si_ref, si)
        end
    end
end
