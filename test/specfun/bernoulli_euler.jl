# SPDX-License-Identifier: MIT

@testset "bernoulli.bernoa" begin
    test_n = [
        2:10...,
        rand(10:100, 2)...,
    ]
    for n = test_n
        @testset "bernoa($n)" begin
            bn_ref = _bernoa(n)
            bn = Specfun.bernoa(n)
            @test isapprox(bn_ref, bn)
        end
    end
end

@testset "bernoulli.bernob" begin
    test_n = [
        2:10...,
        rand(10:100, 2)...,
    ]
    for n = test_n
        @testset "bernob($n)" begin
            bn_ref = _bernob(n)
            bn = Specfun.bernob(n)
            @test isapprox(bn_ref, bn)
        end
    end
end


@testset "euler.eulera" begin
    test_n = [
        2:10...,
        rand(10:100, 2)...,
    ]
    for n = test_n
        @testset "eulera($n)" begin
            en_ref = _eulera(n)
            en = Specfun.eulera(n)
            @test isapprox(en_ref, en)
        end
    end
end

@testset "euler.eulerb" begin
    test_n = [
        2:10...,
        rand(10:100, 2)...,
    ]
    for n = test_n
        @testset "eulerb($n)" begin
            en_ref = _eulerb(n)
            en = Specfun.eulerb(n)
            @test isapprox(en_ref, en)
        end
    end
end
