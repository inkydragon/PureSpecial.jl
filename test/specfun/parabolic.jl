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
    test_x = Float64[
        0.0, -0.0,
        -10:10...,
        rand(3)...,
        -rand(5)...,
        1e5*rand(5)...,
        -1e5*rand(5)...,
    ]
    test_v = Float64[
        0.0, -0.0,
        -5:5...,
        rand(5)...,
        -rand(5)...,
    ]
    
    for x in test_x,
        v in test_v
        vv = v + copysign(1.0, v)
        na = abs(trunc(Int, vv)) + 1
    
        pref_dv, ref_dp = zeros(na), zeros(na)
        dv, dp = zeros(na), zeros(na)
        @testset "pbdv!($x, $v)" begin
            r_pdf, r_pdd = _pbdv!(pref_dv,ref_dp, x,v)
            pdf, pdd = Specfun.pbdv!(dv,dp, x,v)
    
            @test isapprox(pref_dv, dv; nans=true)
            @test isapprox(ref_dp, dp; nans=true)
            @test isapprox(r_pdf, pdf; nans=true)
            @test isapprox(r_pdd, pdd; nans=true)
        end
    end
end
