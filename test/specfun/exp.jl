# SPDX-License-Identifier: MIT

@testset "_e1xb" begin
    test_x = Float64[
        0:10...,
        rand(10)...,
    ]
    
    for x in test_x
        @testset "e1xb($x)" begin
            @test isapprox(_e1xb(x), Specfun.e1xb(x))
        end
    end
end
