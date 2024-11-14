# SPDX-License-Identifier: MIT

@testset "lpmns" begin
    test_m = Int[
        0:4...,
    ]
    test_n = Int[
        1:4...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        -0.4,
        # `!(if abs(x) < 1.0001)`
        3.14,
        # `!(if abs(x) < 1.0001) && abs(x) > 1.1`
        prevfloat(1.1),
    ]

    for m in test_m,
        n in test_n,
        x in test_x
        if (m) > n
            continue
        end
        @testset "lpmns(m=$m, n=$n, x=$x)" begin
            pm_ref, pd_ref = zeros(n+2), zeros(n+2)
            pm, pd = zeros(n+2), zeros(n+2)
            Specfun.lpmns(m, n, x, pm, pd)
            _lpmns(m, n, x, pm_ref, pd_ref)
            # 
            @test isapprox(pm_ref, pm; nans=true)
            @test isapprox(pd_ref, pd; nans=true)
        end   
    end
end

@testset "lqmns" begin
    test_m = Int[
        0:4...,
    ]
    test_n = Int[
        1:4...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        -0.4,
        # `!(if abs(x) < 1.0001)`
        3.14,
        # `!(if abs(x) < 1.0001) && abs(x) > 1.1`
        prevfloat(1.1),
    ]

    for m in test_m,
        n in test_n,
        x in test_x
        @testset "lqmns(m=$m, n=$n, x=$x)" begin
            qm_ref, qd_ref = zeros(n+1), zeros(n+1)
            qm, qd = zeros(n+1), zeros(n+1)
            _lqmns(m, n, x, qm_ref, qd_ref)
            Specfun.lqmns(m, n, x, qm, qd)
            # 
            @test isapprox(qm_ref, qm; nans=true)
            @test isapprox(qd_ref, qd; nans=true)
        end   
    end
end

@testset "_lpmv0" begin
    test_v = Float64[
        1.1, 3.14,
        1:5...,
    ]
    test_m = Int[
        0:5...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        # `!(x >= -0.35)`
        -0.4,
    ]

    for v in test_v,
        x in test_x,
        m in test_m
        @assert abs(x) <= 1
        @testset "_lpmv0(v=$v, m=$m, x=$x)" begin
            p_ref = _lpmv0(v, m, x)
            p = Specfun.lpmv0(v, m, x)
            @test isapprox(p_ref, p)
        end   
    end
end

@testset "lpmv" begin
    test_v = Float64[
        -1:5...,
        1.1, 3.14,
        1:5...,
    ]
    test_m = Int[
        -4:4...,
    ]
    test_x = Float64[
        -1, 1,
        0.0, -0.0,
        # `!(x >= -0.35)`
        -0.4,
    ]

    for v in test_v,
        x in test_x,
        m in test_m
        @assert abs(x) <= 1
        @testset "lpmv(v=$v, m=$m, x=$x)" begin
            p_ref = _lpmv(v, m, x)
            p = Specfun.lpmv(v, m, x)
            @test isapprox(p_ref, p; nans=true)
        end   
    end
end
