# SPDX-License-Identifier: MIT

@testset "hyp1f1/special" begin
    function hyp1f1_test(chg, a, b, z::ComplexF64)
        chg = complex(chg)
        @testset "cchg($a, $b, $z)" begin
            @test isapprox(chg, Specfun.cchg(a, b, z))
            @test isapprox(chg, _cchg(a, b, z))
        end
    end
    function hyp1f1_test(hg, a, b, x::Float64)
        @testset "chgm($a, $b, $x)" begin
            @test isapprox(hg, Specfun.chgm(a, b, x))
            # NOTE: _chgm donot handle spcial inputs
            # @test isapprox(hg, _chgm(a, b, x))
            @test isapprox(_chgm(a, b, x), Specfun.chgm_kernel(a, b, x))
        end
        z = complex(x)
        chg = complex(hg)
        hyp1f1_test(chg, a, b, z)
    end

    for _ in 1:5
        let a=rand(), b=rand(), x=rand(), z=rand(ComplexF64)
            # HG = Inf  for b = negative int
            hyp1f1_test(1.0e300, a, -1.0, x)
            hyp1f1_test(1.0e300, a, -1.0, z)
            # HG = 1    for M(0,b,x)
            hyp1f1_test(1, 0.0, b, x)
            hyp1f1_test(1, 0.0, b, z)
            # HG = 1    for M(a,b,0)
            hyp1f1_test(1, a, b, 0.0)
            hyp1f1_test(1, a, b, -0.0-0.0im)
            # HG = 1 - x/b  for M(-1,b,x)
            hyp1f1_test(1-x/b, -1.0, b, x)
            hyp1f1_test(1-z/b, -1.0, b, z)
            # HG = e^x      for M(a,a,x)
            hyp1f1_test(exp(x), a, a, x)
            hyp1f1_test(exp(z), a, a, z)
            # HG = (1 + x/b)*e^x    for M(b+1,b,x)
            hyp1f1_test((1 + x/b)*exp(x), b+1, b, x)
            hyp1f1_test((1 + z/b)*exp(z), b+1, b, z)
            # HG = [exp(x)-1]/x     for M(1,2,x)
            hyp1f1_test((exp(x)-1)/x, 1.0, 2.0, x)
            hyp1f1_test((exp(z)-1)/z, 1.0, 2.0, z)
        end
    end
end

function hyp1f1_test(a, b, z::ComplexF64)
    @testset "cchg($a, $b, $z)" begin
        @test isapprox(_cchg(a, b, z), Specfun.cchg(a, b, z))
    end
end
function hyp1f1_test(a, b, x::Float64)
    @testset "chgm($a, $b, $x)" begin
        # NOTE: _chgm donot handle spcial inputs
        # @test isapprox(_chgm(a, b, x), Specfun.chgm(a, b, x))
        ref = _chgm(a, b, x)
        res = Specfun.chgm_kernel(a, b, x)

        if isnan(ref)
            @test isequal(ref, res)
        else
            @test isapprox(ref, res)
        end
    end
    z = complex(x)
    hyp1f1_test(a, b, z)
end

"""
Branch Coverage Test Inputs for `chgm`
"""
const TEST_HYP1F1_ABX = [
    # ( a, b, x )
    
    #= if isinteger(a) && a < 0.0 =#
    (-23.0, 1.0, 1.0),

    #= if x < 0.0 =#
    # T1: x < 0.0
    (1.0,  2.0, -1.0),
    (1.0, -2.0, -1.0),
    # F1: x >= 0.0
    (1.0,  2.0, 0.0),
    (1.0, -2.0, 0.0),
    (1.0,  2.0, 1.0),
    (1.0, -2.0, 1.0),

    #= if a >= 2.0 =#
    # T2: a >= 2.0
    (2.0, 1.0, 1.0),
    (3.1, 1.0, 1.0),
    # F2: a < 2.0
    (1.0, 2.0, 1.0),
    (-1.0, 2.0, 1.0),

    # NOTE: 
    #   if x < 0.0:     a -= b
    #   if a >= 2.0:    a -= Int(a) + 1

    #= if (abs(x) <= 30.0 + abs(b)) || (a < 0.0) =#
    # T3
    # T31: (T || F):  (abs(x) <= 30+abs(b)) && a >= 0.0
    # T311=> x >= 0 && a in [0, 2) && (x < 30+abs(b))
    (0.0,  0.0, 0.0),
    (0.0,  1.0, 10.0),
    (1.9,  1.0, 10.0),
    (1.9,  1.0, 31.0),
    # T312=> x < 0 && a < 2 && (abs(x) <= 30+abs(b)) && (b-a) > 0
    (0.1,  1.0, -10.0),
    (-0.1, -1.0, -10.0),
    (-0.1, -1.0, -31.0),
    # T32: (F || T):  (abs(x) > 30+abs(b)) && a < 0
    # T321=> x >= 0 && ...
    (-1.9,  1.0, 42.0),
    (-1.9, -1.0, 42.0),
    # T322=> x < 0 && (b-a) < 0 && ...
    (2.0,  1.0, -42.0),
    (2.0, -1.0, -42.0),

    # F3: (abs(x) > 30+abs(b)) && a >= 0.0
    #= if x0 >= 0.0  && F3 =#
    # T4: x0 >= 0
    (1.9,  1.0, 42.0),
    (1.9, -1.0, 42.0),
    # F4: x0 < 0
    (1.9,  2.0, -42.0),
    (-2.1, -2.0, -42.0),

    #= if (n == 0) =#
    # Always

    #= if (n == 1) =#
    # T5: nl >= 1 ==> a >= 2.0
    (31.0, 1.0, 1.0),
    (31.0, -1.0, 1.0),

    #= if (a0 >= 2.0) =#
    # Same as T2/F2
]

@testset "hyp1f1/Branch Coverage" begin
    for (a,b,x) in TEST_HYP1F1_ABX
        hyp1f1_test(a, b, x)
    end
end

@testset "hyp1f1/cchg" begin
    broken_b = Float64(0xBAD)
    test_abz = [
        #= if b == 0.0 || b == -abs(b) =#
        (1.0, 0.0, complex(1.0)),
        (1.0, -0.0, complex(1.0)),
        
        #= if a == 0.0 || z == 0.0 =#
        (0.0, 1.0, complex(0.0)),
        
        #= if a == -1.0 =#
        (-1.0, 1.0, complex(1.0)),
        (-1.0, -1000.0, complex(1.0)),
        (-1.0, 1.0, -complex(1000.0)),
        
        #= if a == b =#
        (1.0, 1.0, complex(3.14, 2.7)),
        
        #= if (a - b) == 1.0 =#
        (13.0, 12.0, complex(1.0)),
        (-12.0, -13.0, complex(-1.0)),
        
        #= if a == b =#
        (1.0, 2.0, complex(1.0)),
        
        #= if isinteger(a) && a < 0.0 =#
        # T1: 
        (-23.0, 1.0, complex(1.0)),
        
        # F1: !isinteger(a) || a >= 0.0
        (-3.14, broken_b, complex(1.0)),
        (3.14, broken_b, complex(1.0)),
        (23.0, broken_b, complex(1.0)),
        
        #= if (abs(z) < 20.0 + abs(b)) || (a < 0.0) =#
        (314.1, broken_b, complex(31.0, 27.0)),
    ]

    for (a,b,z) in test_abz
        @testset "cchg($a, $b, $z)" begin
            ref = _cchg(a, b, z)
            res = Specfun.cchg(a, b, z)

            if isnan(ref)
                @test isequal(ref, res)
            else
                @test isapprox(ref, res)
            end
        end
    end
end

@testset "chgul" begin
    test_ab = Tuple{Float64, Float64}[
        # (il1,il2) = (T, T)
        (0.0, 1.0),
        (-10, 20),
        # (il1,il2) = (T, F)
        (0.0, 0.0),
        (0.0, 3.5),
        (-10, -20),
        # (il1,il2) = (F, T)
        (1.0, 1.0),
        (1.0, 2.0),
        (10, 20),
        # (il1,il2) = (F, F)
        (1.0, 0.0),
        (1.1, 1.2),
        (1.0, -2.0),
        (10, -20),
        (20, 10),
    ]
    test_x = Float64[
        1:20...,
        rand(5)...,
        (rand(1:1000)*rand(5))...,
    ]
    for (a, b) in test_ab,
        x in test_x
        @testset "chgul(a=$a; b=$b; x=$x)" begin
            hu_ref, id_ref = _chgul(a, b, x)
            hu, id = Specfun.chgul(a, b, x)
            @test isapprox(hu_ref, hu)
            @test isapprox(id_ref, id)
        end
    end
end

@testset "chgus" begin
    test_ab = Tuple{Float64, Float64}[
        # a == b
        (3.14, 3.14),
        # a < b
        (-1, 0.1),
        (-1, -0.1),
        (1.0, 2.1),
        # a > b
        (1, 0.1),
        (1, -0.1),
        (10, 1.1),
        (10, -3.1),
    ]
    test_x = Float64[
        1:20...,
        rand(5)...,
        (rand(1:1000)*rand(5))...,
    ]
    for (a, b) in test_ab,
        x in test_x
        @testset "chgus(a=$a; b=$b; x=$x)" begin
            hu_ref, id_ref = _chgus(a, b, x)
            hu, id = Specfun.chgus(a, b, x)
            @test isapprox(hu_ref, hu)
            @test isapprox(id_ref, id)
        end
    end
end

@testset "chgubi" begin
    test_ab = Tuple{Float64, Float64}[
        # a == b
        (1, 1),
        (7, 7),
        # a < b
        (2, 7),
        (3.14, 7),
        # a > b
        (7, 1),
        (3.14, 1),
        (7, -1),
        # (100, 50),
    ]
    test_x = Float64[
        1:20...,
        rand(5)...,
        (rand(1:1000)*rand(5))...,
    ]
    for (a, b) in test_ab,
        x in test_x
        @testset "chgubi(a=$a; b=$b; x=$x)" begin
            hu_ref, id_ref = _chgubi(a, b, x)
            hu, id = Specfun.chgubi(a, b, x)
            @test isapprox(hu_ref, hu)
            @test isapprox(id_ref, id)
        end
    end
end

@testset "chguit" begin
    test_ab = Tuple{Float64, Float64}[
        # a == b
        (1, 1),
        (7, 7),
        # a < b
        (2, 7),
        (3.14, 7),
        # a > b
        (7, 1),
        (3.14, 1),
        (7, -1),
        # (100, 50),
    ]
    test_x = Float64[
        1:20...,
        rand(5)...,
        (rand(1:1000)*rand(5))...,
    ]
    for (a, b) in test_ab,
        x in test_x
        @testset "chguit(a=$a; b=$b; x=$x)" begin
            hu_ref, id_ref = _chguit(a, b, x)
            hu, id = Specfun.chguit(a, b, x)
            @test isapprox(hu_ref, hu)
            @test isapprox(id_ref, id)
        end
    end
end
