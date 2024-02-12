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

@testset "hyp1f1/chgm" begin

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
