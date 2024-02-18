# SPDX-License-Identifier: MIT

const _GAM0_TEST_X = Float64[
    -1.0:0.01:1.0...,
    -1*rand(10)...,
    0.0, -0.0,
    rand(10)...,
]

@testset "_gam0" begin
    for x in _GAM0_TEST_X
        @testset "gam0($x)" begin
            if -1 == x
                @test_broken false
                continue
            end

            # NOTE: _gam0(x) doesn't give right answer for big x
            if abs(x) >= 0.4
                ref = Specfun.cgama(x, 1)
                @test isapprox(real(ref), Specfun.gam0(x))
            else 
                @test isapprox(_gam0(x), Specfun.gam0(x))
            end
        end   
    end
end

const _GAMMA2_TEST_X = Float64[
    _GAM0_TEST_X...,

    -42:42...,
    rand(-1000:1000, 10)...,
]

@testset "_gamma2" begin
    for x in _GAMMA2_TEST_X
        @testset "gamma2($x)" begin
            @test isapprox(_gamma2(x), Specfun.gamma2(x))
            
            if isinteger(x)
                xp = nextfloat(x)
                @test isapprox(_gamma2(xp), Specfun.gamma2(xp))
                xm = prevfloat(x)
                @test isapprox(_gamma2(xm), Specfun.gamma2(xm))
            end
        end   
    end
end


@testset "_gaih" begin
    test_x = Float64[
        0.1,
        1:100...,
    ]
    test_x /= 2

    for x in test_x
        @testset "gaih($x)" begin
            @test isequal(_gaih(x), Specfun.gaih(x))
        end   
    end
end
