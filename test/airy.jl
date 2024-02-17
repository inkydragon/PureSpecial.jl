# SPDX-License-Identifier: MIT


"""Zeros of `Ai(x)` and `Ai'(x)`
Ref: https://dlmf.nist.gov/9.9#T1
"""
const DLMF_TABLE_9_9_1 = [
#   k,    a_k,           Ai'(a_k),       a'_k,           Ai(a'_k)
    1    -2.33810_74105  0.70121_08227   -1.01879_29716  0.53565_66560;
    2    -4.08794_94441 -0.80311_13697   -3.24819_75822 -0.41901_54780;
    3    -5.52055_98281  0.86520_40259   -4.82009_92112  0.38040_64686;
    4    -6.78670_80901 -0.91085_07370   -6.16330_73556 -0.35790_79437;
    5    -7.94413_35871  0.94733_57094   -7.37217_72550  0.34230_12444;
    6    -9.02265_08533 -0.97792_28086   -8.48848_67340 -0.33047_62291;
    7   -10.04017_43416  1.00437_01227   -9.53544_90524  0.32102_22882;
    8   -11.00852_43037 -1.02773_86888  -10.52766_03970 -0.31318_53910;
    9   -11.93601_55632  1.04872_06486  -11.47505_66335  0.30651_72939;
   10   -12.82877_67529 -1.06779_38592  -12.38478_83718 -0.30073_08293;
]

@testset "ai_zeros" begin
    @test_throws DomainError ai_zeros(-1)

    @test ai_zeros(0) == (Float64[], Float64[], Float64[], Float64[])

    a, da, ai, dai = ai_zeros(10)
    @test isapprox(DLMF_TABLE_9_9_1[:,2] , a)
    @test isapprox(DLMF_TABLE_9_9_1[:,4] , da)
    @test isapprox(DLMF_TABLE_9_9_1[:,5] , ai)
    @test isapprox(DLMF_TABLE_9_9_1[:,3] , dai)
end


"""Zeros of `Bi(x)` and `Bi'(x)`
Ref: https://dlmf.nist.gov/9.9#T2
"""
const DLMF_TABLE_9_9_2 = [
#    k,   b_k,           Bi'(b_k),        b'_k,          Bi(b'_k)
     1   -1.17371_32227  0.60195_78880   -2.29443_96826 -0.45494_43836;
     2   -3.27109_33028 -0.76031_01415   -4.07315_50891  0.39652_28361;
     3   -4.83073_78417  0.83699_10126   -5.51239_57297 -0.36796_91615;
     4   -6.16985_21283 -0.88947_99014   -6.78129_44460  0.34949_91168;
     5   -7.37676_20794  0.92998_36386   -7.94017_86892 -0.33602_62401;
     6   -8.49194_88465 -0.96323_44302   -9.01958_33588  0.32550_97364;
     7   -9.53819_43793  0.99158_63705  -10.03769_63349 -0.31693_46537;
     8  -10.52991_35067 -1.01638_96592  -11.00646_26677  0.30972_59408;
     9  -11.47695_35513  1.03849_42860  -11.93426_16450 -0.30352_76648;
    10  -12.38641_71386 -1.05847_18444  -12.82725_83092  0.29810_49111;
]

@testset "bi_zeros" begin
    @test_throws DomainError bi_zeros(-1)

    @test bi_zeros(0) == (Float64[], Float64[], Float64[], Float64[])
    
    b, db, bi, dbi = bi_zeros(10)
    @test isapprox(DLMF_TABLE_9_9_2[:,2] , b)
    @test isapprox(DLMF_TABLE_9_9_2[:,4] , db)
    @test isapprox(DLMF_TABLE_9_9_2[:,5] , bi)
    @test isapprox(DLMF_TABLE_9_9_2[:,3] , dbi)
end

"""
ref: CoSF Table 10.2: Integrals of Airy Functions

Table Generation:
```wolfram
data=N[
    Table[{x, 
        Integrate[AiryAi[t], {t, 0, x}],  Integrate[AiryBi[t], {t, 0,x}], 
        Integrate[AiryAi[-t], {t, 0, x}],  Integrate[AiryBi[-t], {t, 0,x}]
    }, {x, 0, 10, 0.2}],
18];
fmtNum[n_]:=NumberForm[n, 8]
formattedData=Map[fmtNum,data]
TableForm[formattedData]
```

Note: results given by wolfram slightly differ from the values in the book.
"""
const COSF_TABLE_10_2 = Float32[
#    x,   apt,        bpt,        ant,        bnt
     0.0  0.0         0.0         0.0         0.0;
     0.2  0.065851514 0.13199449  0.076156954 0.11398096;
     0.4  0.12164061  0.28256702  0.16229441  0.20952889;
     0.6  0.16801789  0.45356503  0.25736069  0.28553622;
     0.8  0.20589445  0.64845817  0.35944153  0.3405258;
     1.0  0.23631734  0.87276912  0.46567398  0.37300501;
     1.2  0.26037122  1.1346638   0.57224048  0.38185428;
     1.4  0.27910665  1.4457942   0.67447310  0.36673336;
     1.6  0.29349241  1.8225233   0.76709256  0.3284724;
     1.8  0.30438816  2.2877211   0.84459409  0.26939973;
     2.0  0.31253276  2.8734083   0.90177283  0.19354741;
     2.2  0.31854425  3.6246635   0.93435563  0.10667181;
     2.4  0.32292739  4.6054276   0.93967668  0.016034543;
     2.6  0.32608569  5.9071710   0.91730538 -0.070090123;
     2.8  0.32833548  7.6619246   0.86951369 -0.14317875;
     3.0  0.32992038  10.062003   0.80146284 -0.19544249;
     3.2  0.33102486  13.390077   0.72100368 -0.22092493;
     3.4  0.33178649  18.065367   0.63802561 -0.2165557;
     3.6  0.33230632  24.715120   0.56335609 -0.18296472;
     3.8  0.33265759  34.286066   0.50730045 -0.12485427;
     4.0  0.33289265  48.219475   0.47800750 -0.05076006;
     4.2  0.33304843  68.728194   0.47992943  0.027887893;
     4.4  0.33315072  99.238280   0.51269276  0.098370159;
     4.6  0.33321726  145.09820   0.57068590  0.14876495;
     4.8  0.33326017  214.72601   0.64358513  0.17018589;

     5.0  0.33328759  321.47832   0.71788220  0.15873094;
     5.2  0.33330496  486.71706   0.77926271  0.11667298;
     5.4  0.33331588  744.87753   0.81545487  0.05250028;
     5.6  0.33332268  1151.9023   0.81897895 -0.020389929;
     5.8  0.33332688  1799.3751   0.78914055 -0.086251832;
     6.0  0.33332945  2838.3697   0.73267526 -0.13038106;
     6.2  0.33333102  4519.9696   0.66268959 -0.14262054;
     6.4  0.33333196  7264.5788   0.59592616 -0.12011152;
     6.6  0.33333253  11781.281   0.54883591 -0.068472877;
     6.8  0.33333286  19274.766   0.53334743 -0.00088803208;
     7.0  0.33333306  31806.438   0.55345167  0.064916739;
     7.2  0.33333318  52928.749   0.60365964  0.11121471;
     7.4  0.33333324  88806.742   0.66999959  0.12521801;
     7.6  0.33333328  150213.86   0.73355341  0.10300571;
     7.8  0.33333330  256106.31   0.77575129  0.051143486;
     8.0  0.33333332  440065.25   0.78398260 -0.014756444;
     8.2  0.33333332  761981.09   0.75578555 -0.074404324;
     8.4  0.33333333  1.3293781e6 0.70011697 -0.10902216;
     8.6  0.33333333  2.3365751e6 0.63499076 -0.10749347;
     8.8  0.33333333  4.1370433e6 0.58192697 -0.070396425;
     9.0  0.33333333  7.3779210e6 0.55881975 -0.010330376;
     9.2  0.33333334  1.3251594e7 0.57358506  0.051922358;
     9.4  0.33333334  2.3969060e7 0.62093760  0.09439868;
     9.6  0.33333333  4.3655899e7 0.68375251  0.10183701;
     9.8  0.33333334  8.0058205e7 0.73889845  0.071573317;
    10.0  0.33333332  1.4780980e8 0.76569840  0.015040428;
]

@testset "itairy" begin
    @test_throws DomainError bi_zeros(-1)

    @test isequal(itairy(0.0), (0.0, 0.0, 0.0, 0.0))

    for (x, r_apt, r_bpt, r_ant, r_bnt) in eachrow(COSF_TABLE_10_2)
        @testset "itairy($x)" begin
            apt, bpt, ant, bnt = itairy(Float64(x))

            @test isapprox(r_apt, Float32(apt))
            @test isapprox(r_bpt, Float32(bpt))
            @test isapprox(r_ant, Float32(ant))
            @test isapprox(r_bnt, Float32(bnt))
        end
    end
end
