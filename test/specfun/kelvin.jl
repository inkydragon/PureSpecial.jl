# SPDX-License-Identifier: MIT

@testset "Kelvin functions" begin
    special_inputs = Float64[
    #= if xa == 0.0 =#
        -0.0, 0.0,
        1:20...,
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
        end
    end
end
