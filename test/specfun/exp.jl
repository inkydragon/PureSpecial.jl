# SPDX-License-Identifier: MIT

const E1XB_TEST_X = Float64[
    0:10...,
    rand(10)...,
]

@testset "_e1xb" begin
    for x in E1XB_TEST_X
        @testset "e1xb($x)" begin
            @test isapprox(_e1xb(x), Specfun.e1xb(x))
        end
    end
end

@testset "_e1z" begin
    test_z = ComplexF64[
        E1XB_TEST_X...,
        -50:1:50...,
    ]
    test_z = ComplexF64[
        test_z...,
        (x*im for x in test_z)...,
        (x + rand(-100:100) * im for x in test_z)...,
    ]

    for z in test_z
        @testset "e1z($z)" begin
            @test isapprox(_e1z(z), Specfun.e1z(z))
        end
    end
end
