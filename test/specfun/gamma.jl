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

@testset "lgama" begin
    test_x = Float64[
        1:10...,
        rand(10)...,
    ]
    for x in test_x,
        kf in 0:1
        @testset "lgama(kf=$kf; x=$x)" begin
            ref = _lgama(kf, x)
            res = Specfun.lgama(kf, x)
            @test isapprox(ref, res)
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

@testset "_cgama" begin
    test_z = [
        1.0 + im,
        1.0 - im,
        -1.0 + im,
        -1.0 - im,
        
        -1.0 + 0im,
        -2.0 + 1im,
        2.0 + 3im,
        
        # f77 acc
        0.75 + 2.166936862967124im
    ]

    for z in test_z,
        kf in 0:1
        @testset "cgama($z, $kf)" begin
            ref = _cgama(z, kf)
            res = Specfun.cgama(z, kf)
            @test isapprox(ref, res)
        end
    end
end

@testset "beta" begin
    test_p = Float64[
        1:5...,
        rand(5)...,
    ]
    test_q = Float64[
        1:5...,
        rand(5)...,
    ]
    for p in test_p,
        q in test_q
        @testset "beta(p=$p; q=$q)" begin
            ref = _beta(p, q)
            res = Specfun.beta(p, q)
            @test isapprox(ref, res)
        end
    end
end


@testset "incog" begin
    test_a = Float64[
        1:5...,
        rand(5)...,
        rand(6:170)...,
        200
    ]
    test_x = Float64[
        0:5...,
        rand(5)...,
        rand(6:700)...,
        800
    ]
    for a in test_a,
        x in test_x
        @testset "incog(a=$a; x=$x)" begin
            gin_r, gim_r, gip_r, isfer_r = _incog(a, x)
            gin, gim, gip, isfer = Specfun.incog(a, x)
            @test isapprox(isfer_r, isfer)
            @test isapprox(gin_r, gin; nans=true)
            @test isapprox(gim_r, gim; nans=true)
            @test isapprox(gip_r, gip; nans=true)
        end
    end
end


@testset "psi" begin
    test_x = Float64[
        # negative int: `(x == trunc(Int, x)) && (x <= 0.0)`
        -4:-1...,
        # nonnegative int:  `if xa == trunc(Int, xa)`
        0:10...,
        -0.0,
        # N+0.5:    `(xa + 0.5) == trunc(Int, xa + 0.5)`
        1.5, 2.5, 3.5,
        # `else`
        3.1,
        eps(),
        rand(4)...,
        # `x < 0.0`
        -rand(4)...,
    ]

    for x in test_x
        @testset "psi($x)" begin
            p_ref = _psi(x)
            p = Specfun.psi(x)
            @test isapprox(p_ref, p)
        end   
    end
end

@testset "cpsi" begin
    test_z = [
        0.0 + im,
        1.0 + im,
        1.0 - im,
        -1.0 + im,
        -1.0 - im,
        -1.0 + 0im,
        -2.0 + 1im,
        2.0 + 3im,
    ]
    for z in test_z
        @testset "cpsi($z)" begin
            ref = _cpsi(z)
            res = Specfun.cpsi(z)
            @test isapprox(ref, res)
        end
    end
end
