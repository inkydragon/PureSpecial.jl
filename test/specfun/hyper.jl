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

end
