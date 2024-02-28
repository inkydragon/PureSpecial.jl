# SPDX-License-Identifier: MIT

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
