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
    test_a = Float64[
        -0.0,
        -5:5...,
        rand(-5.0:eps():5.0, 10)...,
        # match f77 acc
        4.333873725934248,
    ]
    test_x = Float64[
        -0.0,
        -5:5...,
        rand(-5.0:eps():5.0, 10)...,
    ]

    for a in test_a,
        x in test_x
        @testset "pbwa($a, $x)" begin
            r_w1f, r_w1d, r_w2f, r_w2d = _pbwa(a, x)
            w1f, w1d, w2f, w2d = Specfun.pbwa(a, x)
            # TODO: fix rtol
            @test isapprox(r_w1f, w1f; nans=true, rtol=1e-7)
            @test isapprox(r_w1d, w1d; nans=true, rtol=1e-7)
            @test isapprox(r_w2f, w2f; nans=true, rtol=1e-7)
            @test isapprox(r_w2d, w2d; nans=true, rtol=1e-7)
        end
    end
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

@testset "pbdv/pbvv" begin
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
    
        r_dv, r_dp = zeros(na), zeros(na)
        dv, dp = zeros(na), zeros(na)
        @testset "pbdv($x, $v)" begin
            r_pdf, r_pdd = _pbdv(r_dv,r_dp, x,v)
            pdf, pdd = Specfun.pbdv(dv,dp, x,v)
    
            @test isapprox(r_dv, dv; nans=true)
            @test isapprox(r_dp, dp; nans=true)
            @test isapprox(r_pdf, pdf; nans=true)
            @test isapprox(r_pdd, pdd; nans=true)
        end
    
        arr_len = max(na+1, 3)
        r_dv, r_dp = zeros(arr_len), zeros(arr_len)
        dv, dp = zeros(arr_len), zeros(arr_len)
        @testset "pbvv($x, $v)" begin
            r_pvf, r_pvd = _pbvv(r_dv,r_dp, x,v)
            pvf, pvd = Specfun.pbvv(dv,dp, x,v)
            
            # TODO: fix rtol
            @test isapprox(r_dv, dv; nans=true, rtol=1e-7)
            @test isapprox(r_dp, dp; nans=true, rtol=1e-7)
            @test isapprox(r_pvf, pvf; nans=true, rtol=1e-7)
            @test isapprox(r_pvd, pvd; nans=true, rtol=1e-7)
        end
    end
end
