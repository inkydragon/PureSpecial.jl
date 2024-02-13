# SPDX-License-Identifier: MIT

@testset "_cgama" begin
    test_z = [
        1.0 + im,
        1.0 - im,
        -1.0 + im,
        -1.0 - im,
        
        -1.0 + 0im,
        -2.0 + 1im,
        2.0 + 3im,
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
            
            if broken_b == b
                @test_broken false
                continue
            end

            @test isapprox(ref, res)
        end
    end
end

@testset "hyp1f1/chgm" begin
    test_abx = [
        # ( a, b, x )
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

        #= if (abs(x) <= 30.0 + abs(b)) || (a < 0.0) =#
        # T31: (T || F):  (abs(x) > 30+abs(b)) && a >= 0.0
        (0.0,  1.0, 10.0),
        (1.9,  1.0, 10.0),
        (1.9, -1.0, 10.0),
        (0.0,  1.0, -10.0),
        (1.9,  1.0, -10.0),
        (1.9, -1.0, -10.0),
        # T32: (F || T):  (abs(x) > 30+abs(b)) && a < 0.0
        (-1.9,  1.0, 42.0),
        (-1.9, -1.0, 42.0),
        (-1.9,  1.0, -42.0),
        (-1.9, -1.0, -42.0),
    
        # F3: (abs(x) > 30+abs(b)) && a >= 0.0
        (1.9,  1.0, 42.0),
        (1.9, -1.0, 42.0),
        (1.9,  1.0, -42.0), # x = abs(x)
        (1.9, -1.0, -42.0),

        #= if x0 >= 0.0  && F3 =#
        # T4: (x > 30+abs(b)) && a >= 0.0
        (1.9,  1.0, 314.0),
        (1.9, -1.0, 314.0),

        # F4: (abs(x) > 30+abs(b)) && a >= 0.0 && x < 0
        (1.9,  1.0, -314.0),
        (1.9, -1.0, -314.0),
    ]

    for (a,b,x) in test_abx
        @testset "chgm($a, $b, $x)" begin
            ref = _chgm(a, b, x)
            res = Specfun.chgm(a, b, x)

            if isnan(ref)
                @test isequal(ref, res)
            else
                @test isapprox(ref, res)
            end
        end
    end
end
