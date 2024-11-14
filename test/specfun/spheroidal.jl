# SPDX-License-Identifier: MIT

"""Spheroidal Wave Functions
- specfun_segv
- specfun_rswfp
- specfun_aswfa
- specfun_rswfo
"""

@testset "segv" begin
    # TODO: test branch: `if k != 1 && h[k] < h[k-1]`
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 20),
        (50, 100),
        (100, 200),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9, eps(), 1e-10,
        1:10...,
        rand(10)...,
    ]

    for (m, n) in test_mn,
        c in test_c,
        kd in [1, -1]
        #
        max_cv_len = n - m + 1
        @assert 0 <= max_cv_len <= 200
        icm = (n - m + 2) ÷ 2
        max_eg_nz_len = icm * 2
        eg_ref = zeros(Float64, 200)
        eg_res = zeros(Float64, 200)
        #
        cv_ref, _ = _segv(m, n, c, kd, eg_ref)
        cv_res, _ = Specfun.segv(m, n, c, kd, eg_res)
        @testset "_segv(m=$m, n=$n, c=$c, kd=$kd)" begin
            # Result
            @test isapprox(cv_ref, cv_res)
            @test isapprox(eg_ref[1:max_cv_len], eg_res[1:max_cv_len])
            # temp value
            @test isapprox(eg_ref[max_cv_len+1:max_eg_nz_len], eg_res[max_cv_len+1:max_eg_nz_len])
            # zeros
            @test iszero(eg_ref[max_eg_nz_len+1:end])
            @test iszero(eg_res[max_eg_nz_len+1:end])
        end
    end
end

@testset "sdmn" begin
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 20),
        (50, 100),
        (100, 200),
        (710, 1000),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9, eps(), 1e-10,
        1:10...,
        rand(10)...,
    ]
    test_cv = Float64[
        rand(10)...,
        1:10...,
        # test branch: `kb > 2 && if abs(f) > T(1e100)`
        930:950...,
    ]

    for (m, n) in test_mn,
        kd in [1, -1],
        c in test_c,
        cv in test_cv
        #
        max_df_len = 1 + 25 + trunc(Int, 0.5*(n-m) + c)
        @assert 0 <= max_df_len <= 200
        df_ref = zeros(Float64, 200)
        df_res = zeros(Float64, 200)
        #
        _sdmn(m, n, c, cv, kd, df_ref)
        Specfun.sdmn!(m, n, c, cv, kd, df_res)
        @testset "_sdmn(m=$m,n=$n, c=$c,cv=$cv, kd=$kd)" begin
            # Result
            @test isapprox(df_ref[1:max_df_len], df_res[1:max_df_len]; nans=true)
            # zeros
            # @test iszero(df_ref[max_df_len+1:end])
            # @test iszero(df_res[max_df_len+1:end])
        end
    end
end

@testset "sckb" begin
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 20),
        (100, 200),
        (710, 1000),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9, eps(), 1e-10,
        1:10...,
        rand(10)...,
    ]
    test_cv = Float64[
        rand(10)...,
        1:10...,
    ]

    for (m, n) in test_mn,
        kd in [1, -1],
        c in test_c,
        cv in test_cv
        #
        df_ref = zeros(Float64, 200)
        ck_ref = zeros(Float64, 200)
        df_res = zeros(Float64, 200)
        ck_res = zeros(Float64, 200)
        #
        _sdmn(m, n, c, cv, kd, df_ref)
        _sckb(m, n, c, df_ref, ck_ref)
        Specfun.sdmn!(m, n, c, cv, kd, df_res)
        Specfun.sckb!(m, n, c, df_res, ck_res)
        @testset "_sckb(m=$m, n=$n, c=$c)" begin
            # Result
            @test isapprox(ck_ref, ck_res; nans=true)
        end
    end
end

@testset "aswfa" begin
    test_mn = Tuple{Int64,Int64}[
        # branch: x == 1
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (10, 20),
        (100, 200),
        (710, 1000),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9,
        rand(5)...,
    ]
    test_cv = Float64[
        rand(5)...,
    ]
    test_x = Float64[
        1, -1,
        rand(4)...,
        -rand(4)...,
    ]

    for (m, n) in test_mn,
        kd in [1, -1],
        c in test_c,
        cv in test_cv,
        x in test_x
        #
        s1f_ref, s1d_ref = _aswfa(m, n, c, x, kd, cv)
        s1f_res, s1d_res = Specfun.aswfa(m, n, c, x, kd, cv)
        @testset "_aswfa(m=$m,n=$n, c=$c,x=$x,kd=$kd,cv=$cv)" begin
            # Result
            @test isapprox(s1f_ref, s1f_res; nans=true)
            @test isapprox(s1d_ref, s1d_res; nans=true)
        end
    end
end

@testset "sphj" begin
    test_n = Int[
        0:10...,
    ]
    test_x = Float64[
        eps(0.0),
        1, -1,
        rand(4)...,
        -rand(4)...,
    ]

    for n in test_n,
        x in test_x
        sj_ref = zeros(Float64, 200)
        dj_ref = zeros(Float64, 200)
        sj_res = zeros(Float64, 200)
        dj_res = zeros(Float64, 200)
        #
        _, _, nm_ref = _sphj(n, x, sj_ref, dj_ref)
        _, _, nm_res = Specfun.sphj!(n, x, sj_res, dj_res)
        @testset "sphj!(n=$n, x=$x)" begin
            # Result
            @test isapprox(nm_ref, nm_res)
            @test isapprox(sj_ref, sj_res; nans=true)
            @test isapprox(dj_ref, dj_res; nans=true)
        end
    end
end

@testset "rmn1" begin
    test_mn = Tuple{Int64,Int64}[
        # branch: x == 1
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (10, 20),
        (100, 200),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        1e-9,
        rand(5)...,
    ]
    test_cv = Float64[
        rand(5)...,
    ]
    test_x = Float64[
        0.0,
        1, -1,
        rand(4)...,
        -rand(4)...,
    ]

    for (m, n) in test_mn,
        c in test_c,
        cv in test_cv,
        x in test_x,
        kd in [1, -1]
        if kd / (x * x) >= 1
            continue  # @assert
        end

        df_ref = zeros(Float64, 200)
        df_res = zeros(Float64, 200)
        #
        _sdmn(m, n, c, cv, kd, df_ref)
        r1f_ref, r1d_ref = _rmn1(m, n, c, x, kd, df_ref)
        Specfun.sdmn!(m, n, c, cv, kd, df_res)
        r1f, r1d = Specfun.rmn1(m, n, c, x, kd, df_res)
        @testset "rmn1(m=$m,n=$n,c=$c,x=$x,kd=$kd)+sdmn(cv=$cv)" begin
            # Result
            @test isapprox(r1f_ref, r1f; nans=true)
            @test isapprox(r1d_ref, r1d; nans=true)
        end
    end
end

@testset "sphy.special" begin
    function sphy_nm(n, x)
        sy = zeros(Float64, 201)
        dy = zeros(Float64, 201)
        _, _, nm = Specfun.sphy!(n, x, sy, dy)
        nm
    end

    test_nm = @NamedTuple{x::Float64, nm::Int}.([
        (1e-16, 16),
        (1e-8, 31),
        (1e-4, 52),
        (0.1, 103),
        (1.0, 147),
        (3.14, 183),
        (4.0, 193),
        (4.6, 199),
        (1000, 200),
    ])
    
    for tup in test_nm
        @testset "sphy_nm(nm=$(tup.nm), x=$(tup.x))" begin
            for n in 0:200
                if n > tup.nm
                    @test sphy_nm(n, tup.x) == tup.nm
                else
                    @test sphy_nm(n, tup.x) == n
                end
            end
        end
    end

end

@testset "sphy" begin
    test_x = Float64[
        0.0, 1,
        rand(10)...,
        rand(0:1000, 10)...,
    ]
    test_n = Int[
        0:20...,
    ]

    for x in test_x,
        n in test_n
        sy_ref = zeros(Float64, 200)
        dy_ref = zeros(Float64, 200)
        sy = zeros(Float64, 200)
        dy = zeros(Float64, 200)
        nz_len = n + 1
        #
        _, _, nm_ref = _sphy(n, x, sy_ref, dy_ref)
        _, _, nm = Specfun.sphy!(n, x, sy, dy)
        @testset "_sphy(n=$n, x=$x)" begin
            @test nm_ref == nm
            # Result
            @test isapprox(sy_ref[1:nz_len], sy[1:nz_len]; nans=true)
            @test isapprox(dy_ref[1:nz_len], dy[1:nz_len]; nans=true)
            # zeros
            @test iszero(sy[nz_len+1:end])
            @test iszero(dy[nz_len+1:end])
        end
    end
end

@testset "sdmn" begin
    test_mn = Tuple{Int64,Int64}[
        (1, 2),
        (10, 11),
        (50, 100),
    ]
    test_c = Float64[
        # test br: c < 1e-10
        # 1e-9, eps(), 1e-10,
        1:10...,
        # rand(10)..., # TODO:err
    ]
    test_cv = Float64[
        4.5,
        rand(10)...,
        1:10...,
        # test branch: `kb > 2 && if abs(f) > T(1e100)`
        930:950...,
    ]
    test_x = Float64[
        3.14,
        2:6...,
    ]

    for (m, n) in test_mn,
        c in test_c,
        cv in test_cv,
        x in test_x
        kd = 1
        #
        max_df_len = 1 + 25 + trunc(Int, 0.5*(n-m) + c)
        @assert 0 <= max_df_len <= 200
        @testset "_rmn2l(m=$m,n=$n, c=$c,x=$x, kd=$kd)" begin
            df_ref = zeros(Float64, 200)
            df = zeros(Float64, 200)
            #
            _sdmn(m, n, c, cv, kd, df_ref)
            r2f_ref, r2d_ref, id_ref = _rmn2l(m, n, c, x, kd, df_ref)
            Specfun.sdmn!(m, n, c, cv, kd, df)
            r2f, r2d, id = Specfun.rmn2l(m, n, c, x, kd, df)
            # Result
            @test isapprox(r2f_ref, r2f; nans=true)
            @test isapprox(r2d_ref, r2d; nans=true)
            @test id_ref == id
        end
    end
end
