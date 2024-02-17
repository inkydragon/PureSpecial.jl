# SPDX-License-Identifier: MIT

@testset "_gam0" begin
    test_x = Float64[
        -1.0:0.01:1.0...,
        -1*rand(10)...,
        0.0, -0.0,
        rand(10)...,
    ]
    
    for x in test_x
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
