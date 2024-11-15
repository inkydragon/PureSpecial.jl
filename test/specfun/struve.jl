# SPDX-License-Identifier: MIT

@testset "itsh0" begin
    test_x = Float64[
        eps(0.0),
        rand(4)...,
        0:4...,
        nextfloat(30.0),
        30:34...,
    ]

    for x in test_x
        
        th0_ref = _itsh0(x)
        th0 = Specfun.itsh0(x)
        @testset "itsh0(x=$x)" begin
            @test isapprox(th0_ref, th0; nans=true)
        end
    end
end

@testset "itth0" begin
    test_x = Float64[
        eps(0.0),
        rand(4)...,
        0:4...,
        prevfloat(24.5),
        24.5,
        24:28...,
    ]
    for x in test_x
        th0_ref = _itth0(x)
        th0 = Specfun.itth0(x)
        @testset "itth0(x=$x)" begin
            @test isapprox(th0_ref, th0; nans=true)
        end
    end
end
