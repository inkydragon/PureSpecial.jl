# SPDX-License-Identifier: MIT
"""
- _airyb
- airyzo
- itairy
"""

@testset "airyb" begin
    special_inputs = Float64[
    #= if xa == 0.0 =#
        -0.0, 0.0,

    #= if xa <= xm  && xa != 0.0 =#
    #   T1: `xa <= xm`
    #   T11: x > 0, abs(x) <= 5.0
    #   T12: x < 0, abs(x) <= 8.0
    #   [T1]: x in [-8, 0) ∪ (0, 5]
        -8.0, -eps(-1.0), -eps(0.0),
        eps(0.0), eps(1.0), 5.0,
    #   F1: `xa <= xm` ==> `xa > xm`
    #   F11: x > 0, abs(x) > 5.0
    #   F12: x < 0, abs(x) > 8.0
    #   [F1]: x in (, -8) ∪ (5, )
        -9.0, prevfloat(-8.0),
        nextfloat(5.0), 6.0, 8.0,

    #= if xa < 6.0  && [F1] =#
    #   x in ∅ ∪ (5.0, 6.0)
        nextfloat(5.0), 5.5, prevfloat(6.0),

    #= if xa > 15.0  && [F1] =#
    #   x in (, -15.0) ∪ (15.0, )
        -16.0, prevfloat(-15.0),
        nextfloat(15.0), 16.0,

    #= if x > 0.0  && [F1] =#
    #   T2: x > 0.0
    #   [T2]: x in (5, )
        nextfloat(5.0), 5.1,
    #   F2: x <= 0.0
    #   [F2]: x in (, -8)
        prevfloat(-8.0), -8.1, -70.0,

    #= if xa > 70.0  && [F2] =#
    #   x in (, -70)
        prevfloat(-70.0),
        -500.0,
    #= if xa > 500.0  && [F2] =#
    #   x in (, -500)
        prevfloat(-500.0),
        -1000.0,
    #= if xa > 1000.0  && [F2] =#
    #   x in (, -1000)
        prevfloat(-1000.0),
        -1001.0,

        -150.0,
    #= if xa > 150.0  && [F2] =#
    #   x in (, -150)
        prevfloat(-150.0),
        -3000,
    #= if xa > 3000.0  && [F2] =#
    #   x in (, -3000)
        prevfloat(-3000.0),
        -3001.0,
    ] # special_inputs

    for x in special_inputs
        @testset "_airyb($x)" begin
            ref_ai, ref_bi, ref_ad, ref_bd = _airyb(x)
            ai, bi, ad, bd = Specfun.airyb(x)
            @test isapprox(ref_ai, ai)
            @test isapprox(ref_bi, bi)
            @test isapprox(ref_ad, ad)
            @test isapprox(ref_bd, bd)

            broken_list = [
                -8.0, 5.0, 5.1,
                -500.0
            ]
            if x in broken_list
                @test_broken false
                continue
            end

            @test isequal(ref_ai, ai)
            @test isequal(ref_bi, bi)
            @test isequal(ref_ad, ad)
            @test isequal(ref_bd, bd)
        end
    end
end


@testset "airyzo!" begin
    test_nt = [
        0,
        1:10...,
        rand(10:100)...,
    ]

    for nt in test_nt,
        kf in 1:2
        ra,rb,rc,rd = zeros(nt),zeros(nt),zeros(nt),zeros(nt)
        a,b,c,d = zeros(nt),zeros(nt),zeros(nt),zeros(nt)
        @testset "_airyzo(nt=$nt, kf=$kf)" begin
            _airyzo(nt, kf, ra,rb,rc,rd)
            Specfun.airyzo!(nt, kf, a,b,c,d)
            @test isapprox(ra, a)
            @test isapprox(rb, b)
            @test isapprox(rc, c)
            @test isapprox(rd, d)
        end
    end
end

@testset "itairy" begin
    special_inputs = Float64[
        rand(10)...,
        rand(0.0:eps():9.25, 10)...,
        0:42...,
        rand(1:10^6, 10)...,
        # f77 acc
        8.995242820162433,
        9.03670554040598,
    ]

    for x in special_inputs
        @testset "itairy($x)" begin
            ref_apt, ref_bpt, ref_ant, ref_bnt = _itairy(x)
            apt, bpt, ant, bnt = Specfun.itairy(x)

            @test isapprox(ref_apt, apt)
            @test isapprox(ref_bpt, bpt)
            @test isapprox(ref_ant, ant)
            @test isapprox(ref_bnt, bnt)
        end
    end
end
