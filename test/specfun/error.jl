# SPDX-License-Identifier: MIT

const ERF_TEST_X = Float64[
    -5:5...,
    rand(0.0:eps():3.5, 5)...,
]

@testset "erf" begin
    for x in ERF_TEST_X
        @testset "erf($x)" begin
            r_err = _erf(x)
            err = Specfun.erf(x)

            @test isapprox(r_err, err; nans=true)
        end
    end
end

const ERF_TEST_Z = ComplexF64[
    complex.(ERF_TEST_X)...,
    3.0 + 4.0*im,
    -11.0 - 13.0*im,
]

@testset "erf" begin
    for z in ERF_TEST_Z
        @testset "erf($z)" begin
            r_err = _erf(z)
            err = Specfun.erf(z)

            @test isapprox(r_err, err; nans=true)
        end
    end
end

@testset "cerf" begin
    for z in ERF_TEST_Z
        @testset "cerf($z)" begin
            r_err, r_der = _cerf(z)
            err, der = Specfun.cerf(z)

            @test isapprox(r_err, err; nans=true)
            @test isapprox(r_der, der; nans=true)
        end
    end
end

@testset "cerzo!" begin
    for n in 1:10
        ref_zo = zeros(ComplexF64, n)
        zo = zeros(ComplexF64, n)
        @testset "cerzo!($n)" begin
            _cerzo!(ref_zo, n)
            Specfun.cerzo!(zo, n)

            @test isapprox(ref_zo, zo; nans=true)
        end
    end
end


const CFC_TEST_Z = ComplexF64[
    #= if z==0 =#
    0.0, -0.0,
    #= elseif w0 <= 2.5 =#
    1:2...,
    rand(nextfloat(0.0):eps():2.5, 5)...,
    #= elseif w0 < 4.5 =#
    3:4...,
    rand(nextfloat(2.5):eps():4.5, 5)...,
    #= else =#
    5:10...,
    rand(4.5:eps():42.0, 5)...,
    42+0im, -42+0im,
    42im, -42im,
    -42+42im, 42+42im,
    -42-42im, 42-42im,
]

@testset "cfc" begin
    for z in CFC_TEST_Z
        @testset "cfc(z=$z)" begin
            r_zf, r_zd = _cfc(z)
            zf, zd = Specfun.cfc(z)
            
            if z==4.0
                @test_broken false
                continue 
            end
            
            @test isapprox(r_zf, zf; nans=true)
            @test isapprox(r_zd, zd; nans=true)
        end
    end
end

@testset "cfs" begin
    for z in CFC_TEST_Z
        @testset "cfs($z)" begin
            r_zf, r_zd = _cfs(z)
            zf, zd = Specfun.cfs(z)
            
            if z==4.0
                @test_broken false
                continue 
            end
            
            @test isapprox(r_zf, zf; nans=true)
            @test isapprox(r_zd, zd; nans=true)
        end
    end
end

@testset "fcs" begin
    for x in real.(CFC_TEST_Z)
        @testset "fcs($x)" begin
            r_c, r_s = _fcs(x)
            c, s = Specfun.fcs(x)

            @test isapprox(r_c, c; nans=true)
            @test isapprox(r_s, s; nans=true)
        end
    end
end
