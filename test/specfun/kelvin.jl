# SPDX-License-Identifier: MIT

@testset "Kelvin functions" begin
    special_inputs = Float64[
        -0.0, 0.0,
        1:100...,
    ] # special_inputs

    for x in special_inputs
        @testset "_klvna($x)" begin
            r_ber, r_bei, r_ger, r_gei, r_der, r_dei, r_her, r_hei = _klvna(x)
            ber, bei, ger, gei, der, dei, her, hei = Specfun.klvna(x)

            @test isapprox(r_ber, ber)
            @test isapprox(r_bei, bei)
            @test isapprox(r_ger, ger)
            @test isapprox(r_gei, gei)
            @test isapprox(r_der, der)
            @test isapprox(r_dei, dei)
            @test isapprox(r_her, her)
            @test isapprox(r_hei, hei)

            broken_list = [
                4:9...,
                30, 33, 35,
                42, 45,
                57, 59,
                60,
                71,
                80, 82, 86,
                91,
            ]
            if x in broken_list
                @test_broken false
                continue
            end

            @test isequal(r_ber, ber)
            @test isequal(r_bei, bei)
            @test isequal(r_ger, ger)
            @test isequal(r_gei, gei)
            @test isequal(r_der, der)
            @test isequal(r_dei, dei)
            @test isequal(r_her, her)
            @test isequal(r_hei, hei)
        end
    end
end

@testset "klvnzo" begin
    test_nt = Int64[
        0:20...,
    ]

    for nt in test_nt,
        kd in 1:8
        @testset "_klvnzo(nt=$nt, kd=$kd)" begin
            ref = _klvnzo(nt, kd)
            res = Specfun.klvnzo(nt, kd)
            
            @test isapprox(ref, res)
        end
    end
end
