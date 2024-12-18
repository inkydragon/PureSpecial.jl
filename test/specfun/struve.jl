# SPDX-License-Identifier: MIT

@testset "struve.stvh0" begin
    test_x = Float64[
        eps(0.0),
        # typemax(Float64),
        # NaN, 
        rand(4)...,
        0:4...,
        nextfloat(20.0),
        20:24...,
        1e4,
    ]
    for x in test_x
        ret_ref = _stvh0(x)
        ret = Specfun.stvh0(x)
        @testset "stvh0($x)" begin
            @test isapprox(ret_ref, ret; nans=true)
        end
    end
end

@testset "struve.stvh1" begin
    test_x = Float64[
        eps(0.0),
        # typemax(Float64),
        # NaN, 
        rand(4)...,
        0:4...,
        nextfloat(20.0),
        20:24...,
        1e4,
    ]
    for x in test_x
        ret_ref = _stvh1(x)
        ret = Specfun.stvh1(x)
        @testset "stvh1($x)" begin
            @test isapprox(ret_ref, ret; nans=true)
        end
    end
end

@testset "struve.stvhv" begin
    test_v = Float64[
        -4:3...,
        # -8:12..., 12.5,
        # 0.5, 1.5,
        # rand(5)...,
    ]
    test_x = Float64[
        0:20...,
        # eps(0.0),
        # # typemax(Float64),
        # # NaN, 
        # rand(4)...,
        # 0:4...,
        # nextfloat(20.0),
        # 20:24...,  # TODO: test > 20
        # 1e4,
    ]
    for v in test_v,
        x in test_x
        ret_ref = _stvhv(v, x)
        ret = Specfun.stvhv(v, x)
        @testset "stvhv(v=$v, x=$x)" begin
            @test isapprox(ret_ref, ret; nans=true)
        end
    end
end

@testset "struve.stvl0" begin
    test_x = Float64[
        eps(0.0),
        # typemax(Float64),
        # NaN, 
        rand(4)...,
        0:4...,
        nextfloat(20.0),
        20:24...,
        1e4,
    ]
    for x in test_x
        ret_ref = _stvl0(x)
        ret = Specfun.stvl0(x)
        @testset "stvl0($x)" begin
            @test isapprox(ret_ref, ret; nans=true)
        end
    end
end


@testset "struve.itsh0" begin
    test_x = Float64[
        eps(0.0),
        rand(4)...,
        0:4...,
        nextfloat(30.0),
        30:34...,
        1e4,
    ]
    for x in test_x
        th0_ref = _itsh0(x)
        th0 = Specfun.itsh0(x)
        @testset "itsh0(x=$x)" begin
            @test isapprox(th0_ref, th0; nans=true)
        end
    end
end

@testset "struve.itth0" begin
    test_x = Float64[
        eps(0.0),
        rand(4)...,
        0:4...,
        prevfloat(24.5),
        24.5,
        24:28...,
        1e4,
    ]
    for x in test_x
        th0_ref = _itth0(x)
        th0 = Specfun.itth0(x)
        @testset "itth0(x=$x)" begin
            @test isapprox(th0_ref, th0; nans=true)
        end
    end
end

@testset "struve.itsl0" begin
    test_x = Float64[
        eps(0.0),
        rand(4)...,
        0:4...,
        nextfloat(20.0),
        20:24...,
        1e4,
    ]
    for x in test_x
        th0_ref = _itsl0(x)
        th0 = Specfun.itsl0(x)
        @testset "itsl0(x=$x)" begin
            @test isapprox(th0_ref, th0; nans=true)
        end
    end
end
