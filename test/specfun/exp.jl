# SPDX-License-Identifier: MIT

const E1XB_TEST_X = Float64[
    0:10...,
    rand(10)...,
]

@testset "e1xa" begin
    for x in E1XB_TEST_X
        @testset "e1xa($x)" begin
            @test isapprox(_e1xa(x), Specfun.e1xa(x))
        end
    end
end

@testset "_e1xb" begin
    for x in E1XB_TEST_X
        @testset "e1xb($x)" begin
            @test isapprox(_e1xb(x), Specfun.e1xb(x))
        end
    end
end

const E1Z_TEST_Z = let test_z = E1XB_TEST_X
    test_z = ComplexF64[
        test_z...,
        -50:1:50...,
    ]
    
    ComplexF64[
        test_z...,
        (x*im for x in test_z)...,
        (x + rand(-100:100) * im for x in test_z)...,
    ]
end

@testset "_e1z" begin
    for z in E1Z_TEST_Z
        @testset "e1z($z)" begin
            @test isapprox(_e1z(z), Specfun.e1z(z))
        end
    end
end


@testset "eix" begin
    test_x = [
        E1XB_TEST_X...,
        10:40...,
        rand(4:1000, 5)...,
    ]
    test_x = [
        test_x...,
        -test_x...,
    ]
    for x in test_x
        @testset "eix($x)" begin
            @test isapprox(_eix(x), Specfun.eix(x))
        end
    end
end


@testset "eixz" begin
    for z in E1Z_TEST_Z
        @testset "eixz($z)" begin
            @test isapprox(_eixz(z), Specfun.eixz(z))
        end
    end
end
