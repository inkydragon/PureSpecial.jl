# SPDX-License-Identifier: MIT

"""Parabolic Cylinder functions

- [ ] pbwa
    - [x] cgama

- [ ] pbvv
    - [ ] _gamma2
    - [ ] _vvla
    - [ ] _vvsa

- [ ] pbdv
    - [ ] _dvsa
    - [ ] _dvla

"""

@testset "pbwa" begin
    
end

@testset "_vvla/_dvla" begin
    test_x = Float64[
        1e5*rand(5)...,
        -1e5*rand(5)...,
    ]
    test_va = Float64[
        -5:5...,
    ]
    
    for x in test_x,
        va in test_va
        @testset "vvla($x, $va)" begin
            @test isapprox(_vvla(x,va), Specfun.vvla(x,va); nans=true)
        end
        @testset "dvla($x, $va)" begin
            @test isapprox(_dvla(x,va), Specfun.dvla(x,va); nans=true)
        end
    end
end

@testset "pbvv" begin
    
end

@testset "_dvsa/_vvsa" begin
    test_x = Float64[
        0.0, -0.0,
        -10:10...,
        rand(5)...,
        -rand(5)...,
    ]
    test_va = Float64[
        0.0, -0.0,
        -5:5...,
    ]
    
    for x in test_x,
        va in test_va
        @testset "dvsa($x, $va)" begin
            @test isapprox(_dvsa(x,va), Specfun.dvsa(x,va); nans=true)
        end
        @testset "vvsa($x, $va)" begin
            @test isapprox(_vvsa(x,va), Specfun.vvsa(x,va); nans=true)
        end
    end
end

@testset "pbdv" begin
    
end
