# SPDX-License-Identifier: MIT

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

@testset "euler" begin

end
