# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""
implement:
+ airyzo
    + airyb
+ itairy

not-impl:
- AIRYA(X,AI,BI,AD,BD)
    - AJYIK(X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
"""

"""
    airyb(x::Float64)

Compute Airy functions and their derivatives.

## Input
- `x`: Argument of Airy function

## Output
- (ai, bi, ad, bd)
    - `ai`: Ai(x)
    - `bi`: Bi(x)
    - `ad`: Ai'(x)
    - `bd`: Bi'(x)
"""
function airyb(x::Float64)
    eps = 1.0e-15
    @assert Float64(pi) === 3.141592653589793
    c1 = 0.355028053887817
    c2 = 0.258819403792807
    sr3 = 1.732050807568877

    ai = 0.0
    bi = 0.0
    ad = 0.0
    bd = 0.0

    ck = zeros(Float64, 52)
    dk = zeros(Float64, 52)

    xa = abs(x)
    xq = sqrt(xa)
    xm = 8.0
    if x > 0.0
        xm = 5.0
    end
    if x == 0.0
        ai = c1
        bi = sr3 * c1
        ad = -c2
        bd = sr3 * c2
        return ai, bi, ad, bd
    end

    if xa <= xm
        fx = 1.0
        r = 1.0
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k - 1.0) * x
            fx += r
            abs(r) < (abs(fx) * eps) && break
        end #= 10 =#

        #= 15 =#
        gx = x
        r = x
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k + 1.0) * x
            gx += r
            abs(r) < (abs(gx) * eps) && break
        end #= 20 =#
        #= 25 =#
        ai = c1 * fx - c2 * gx
        bi = sr3 * (c1 * fx + c2 * gx)

        df = 0.5 * x * x
        r = df
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k + 2.0) * x
            df += r
            abs(r) < (abs(df) * eps) && break
        end #= 30 =#

        #= 35 =#
        dg = 1.0
        r = 1.0
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k - 2.0) * x
            dg += r
            abs(r) < (abs(dg) * eps) && break
        end #= 40 =#

        #= 45 =#
        ad = c1 * df - c2 * dg
        bd = sr3 * (c1 * df + c2 * dg)
    else
        @assert xa > xm

        kmax = -1
        km2 = -1
        km = trunc(Int, 24.5 - xa)
        if xa < 6.0
            km = 14
        end
        if xa > 15.0
            km = 10
        end

        if x > 0.0
            kmax = km
        else
            #
            # Choose cutoffs so that the remainder term in asymptotic
            # expansion is epsilon size. The X<0 branch needs to be fast
            # in order to make AIRYZO efficient
            #
            km2 = km
            
            if xa > 70.0
                km = 3
            end
            if xa > 500.0
                km = 2
            end
            if xa > 1000.0
                km = 1
            end

            km2 = km
            if xa > 150.0
                km2 = 1
            end
            if xa > 3000.0
                km2 = 0
            end
            kmax = 2 * km + 1
        end # x <=> 0.0

        xe = xa * xq / 1.5
        xr1 = 1.0 / xe
        xar = 1.0 / xq
        xf = sqrt(xar)
        rp = 0.5641895835477563
        r = 1.0
        for k = 1:kmax
            r *= (6.0 * k - 1.0) / 
                216.0 * (6.0 * k - 3.0) / 
                k * (6.0 * k - 5.0) / 
                (2.0 * k - 1.0)
            ck[k] = r
            dk[k] = -(6.0 * k + 1.0) / (6.0 * k - 1.0) * r
        end #= 50 =#

        if x > 0.0
            sai = 1.0
            sad = 1.0
            r = 1.0
            for k = 1:km
                r *= -xr1
                sai += ck[k] * r
                sad += dk[k] * r
            end #= 55 =#
            
            sbi = 1.0
            sbd = 1.0
            r = 1.0
            for k = 1:km
                r *= xr1
                sbi += ck[k] * r
                sbd += dk[k] * r
            end #= 60 =#

            xp1 = exp(-xe)
            ai = 0.5 * rp * xf * xp1 * sai
            bi = rp * xf / xp1 * sbi
            ad = -0.5 * rp / xf * xp1 * sad
            bd = rp / xf / xp1 * sbd
        else
            @assert x <= 0.0

            xcs = cos(xe + pi / 4.0)
            xss = sin(xe + pi / 4.0)

            ssa = 1.0
            sda = 1.0
            r = 1.0
            xr2 = 1.0 / (xe * xe)
            for k = 1:km
                r *= -xr2
                ssa += ck[2*k] * r
                sda += dk[2*k] * r
            end #= 65 =#

            ssb = ck[1] * xr1
            sdb = dk[1] * xr1
            r = xr1
            for k = 1:km2
                r *= -xr2
                ssb += ck[2*k + 1] * r
                sdb += dk[2*k + 1] * r
            end #= 70 =#

            ai =  rp * xf * (xss*ssa - xcs*ssb)
            bi =  rp * xf * (xcs*ssa + xss*ssb)
            ad = -rp / xf * (xcs*sda + xss*sdb)
            bd =  rp / xf * (xss*sda - xcs*sdb)
        end # x <=> 0.0
    end # xa <=> xm

    return ai, bi, ad, bd
end

"""
    airyzo!(
        nt::Int, kf::Int, 
        xa::Vector{Float64}, xb::Vector{Float64}, xc::Vector{Float64}, xd::Vector{Float64}
    )

Compute the first NT zeros of Airy functions
Ai(x) and Ai'(x), a and a', and the associated
values of Ai(a') and Ai'(a); and the first NT
zeros of Airy functions Bi(x) and Bi'(x), b and
b', and the associated values of Bi(b') and
Bi'(b).

## Example
```jl
nt = 4;
a,b,c,d = zeros(nt),zeros(nt),zeros(nt),zeros(nt)
airyzo!(nt, 1, a,b,c,d)
@show a b c d;
```

## Input
- NT --- Total number of zeros
- KF --- Function code
    - KF=1 for Ai(x) and Ai'(x)
    - KF=2 for Bi(x) and Bi'(x)

## Output
- XA(m) --- a, the m-th zero of Ai(x) or
            b, the m-th zero of Bi(x)
- XB(m) --- a', the m-th zero of Ai'(x) or
            b', the m-th zero of Bi'(x)
- XC(m) --- Ai(a') or Bi(b')
- XD(m) --- Ai'(a) or Bi'(b)
            ( m --- Serial number of zeros )

## Routine called
AIRYB for computing Airy functions and their derivatives
"""
function airyzo!(
    nt::Int, kf::Int, 
    xa::Vector{Float64}, xb::Vector{Float64},
    xc::Vector{Float64}, xd::Vector{Float64})
    @assert Float64(pi) === 3.141592653589793
    
    # TODO: check params

    # move local var
    u = 0.0
    u1 = 0.0
    ai = 0.0
    bi = 0.0
    ad = 0.0
    bd = 0.0
    err = 0.0

    rt = 0.0
    for i = 1:nt
        rt0 = 0.0
        if kf == 1
            u = 3.0 * pi * (4.0 * i - 1) / 8.0
            u1 = 1 / (u * u)
        elseif kf == 2
            if i == 1
                rt0 = -1.17371
            else
                u = 3.0 * pi * (4.0 * i - 3.0) / 8.0
                u1 = 1 / (u * u)
            end
        end

        if rt0 == 0
            #
            # DLMF 9.9.18
            #
            rt0 = -(u*u)^(1.0/3.0) * (
                    1.0
                    + u1 * (5.0 / 48.0
                    + u1 * (-5.0 / 36.0
                    + u1 * (77125.0 / 82944.0
                    + u1 * (-108056875.0 / 6967296.0)))))
        end

        while true
            #= 10 =#
            x = rt0
            ai, bi, ad, bd = airyb(x)

            if kf == 1
                rt = rt0 - ai / ad
            elseif kf == 2
                rt = rt0 - bi / bd
            end

            err = abs((rt - rt0) / rt)
            if err <= 1.0e-12
                break
            else
                rt0 = rt
            end
        end

        xa[i] = rt
        if err > 1.0e-14
            ai, bi, ad, bd = airyb(rt)
        end

        if kf == 1
            xd[i] = ad
        elseif kf == 2
            xd[i] = bd
        end
    end #= 15 =## for i = 1:nt

    for i = 1:nt
        rt0 = 0.0
        if kf == 1
            if i == 1
                rt0 = -1.01879
            else
                u = 3.0 * pi * (4.0 * i - 3.0) / 8.0
                u1 = 1 / (u * u)
            end
        elseif kf == 2
            if i == 1
                rt0 = -2.29444
            else
                u = 3.0 * pi * (4.0 * i - 1.0) / 8.0
                u1 = 1 / (u * u)
            end
        end

        if rt0 == 0
            #
            # DLMF 9.9.19
            #
            rt0 = -(u*u)^(1.0/3.0) * (
                    1.0
                    + u1 * (-7.0 / 48.0
                    + u1 * (35.0 / 288.0
                    + u1 * (-181223.0 / 207360.0
                    + u1 * (18683371.0 / 1244160.0)))))
        end

        while true
            #= 20 =#
            x = rt0
            ai, bi, ad, bd = airyb(x)

            if kf == 1
                rt = rt0 - ad / (ai * x)
            elseif kf == 2
                rt = rt0 - bd / (bi * x)
            end

            err = abs((rt - rt0) / rt)
            if err <= 1.0e-12
                break
            else
                rt0 = rt
            end
        end
        xb[i] = rt

        if err > 1.0e-14
            ai, bi, ad, bd = airyb(rt)
        end

        if kf == 1
            xc[i] = ai
        elseif kf == 2
            xc[i] = bi
        end
    end #= 25 =## for i = 1:nt
end

"""
    itairy(x::Float64)

Compute the integrals of Airy fnctions with
respect to t from 0 and x ( x ≥ 0 )

## Input
x   --- Upper limit of the integral

## Output
apt --- Integration of Ai(t) from 0 and x
bpt --- Integration of Bi(t) from 0 and x
ant --- Integration of Ai(-t) from 0 and x
bnt --- Integration of Bi(-t) from 0 and x
"""
function itairy(x::Float64)
    @assert x >= 0

    # Constants
    EPS = 1e-15
    C1 = 0.355028053887817
    C2 = 0.258819403792807
    SR3 = 1.732050807568877

    # 1/3
    Q0 = 0.3333333333333333
    # 2/3
    Q1 = 0.6666666666666667
    # sqrt(2)
    Q2 = 1.414213562373095
    _A = [
        0.569444444444444e+00, 0.891300154320988e+00,
        0.226624344493027e+01, 0.798950124766861e+01,
        0.360688546785343e+02, 0.198670292131169e+03,
        0.129223456582211e+04, 0.969483869669600e+04,
        0.824184704952483e+05, 0.783031092490225e+06,
        0.822210493622814e+07, 0.945557399360556e+08,
        0.118195595640730e+10, 0.159564653040121e+11,
        0.231369166433050e+12, 0.358622522796969e+13,
    ]

    apt, bpt, ant, bnt = 0.0, 0.0, 0.0, 0.0
    if x == 0.0
        return apt, bpt, ant, bnt
    else
        if abs(x) <= 9.25
            ant = 0.0
            bnt = 0.0
            for l = 0:1
                x *= (-1)^l
                fx = x
                r = x
                for k = 1:40
                    r *= (3.0 * k - 2.0) / (3.0 * k + 1.0) * x / (3.0 * k) * x / (3.0 * k - 1.0) * x
                    fx += r
                    if abs(r) < (abs(fx) * EPS)
                        break
                    end
                end
                gx = 0.5 * x * x
                r = gx
                for k = 1:40
                    r *= (3.0 * k - 1.0) / (3.0 * k + 2.0) * x / (3.0 * k) * x / (3.0 * k + 1.0) * x
                    gx += r
                    if abs(r) < (abs(gx) * EPS)
                        break
                    end
                end

                ant = C1*fx - C2*gx
                bnt = SR3 * (C1*fx + C2*gx)
                if l == 0
                    apt = ant
                    bpt = bnt
                else
                    ant = -ant
                    bnt = -bnt
                    x = -x
                end
            end
        else
            xe = x * sqrt(x) / 1.5
            xp6 = 1.0 / sqrt(6.0 * pi * xe)
            su1 = 1.0
            r = 1.0
            xr1 = 1.0 / xe
            for k = 1:16
                r *= -xr1
                su1 += _A[k] * r
            end
            su2 = 1.0
            r = 1.0
            for k = 1:16
                r *= xr1
                su2 += _A[k] * r
            end
            apt = Q0 - exp(-xe) * xp6 * su1
            bpt = 2.0 * exp(xe) * xp6 * su2
            su3 = 1.0
            r = 1.0
            xr2 = 1.0 / (xe * xe)
            for k = 1:8
                r *= -xr2
                su3 += _A[2 * k] * r
            end
            su4 = _A[1] * xr1
            r = xr1
            for k = 1:7
                r *= -xr2
                su4 += _A[2 * k + 1] * r
            end
            su5 = su3 + su4
            su6 = su3 - su4
            ant = Q1 - Q2 * xp6 * (su5 * cos(xe) - su6 * sin(xe))
            bnt = Q2 * xp6 * (su5 * sin(xe) + su6 * cos(xe))
        end
    end

    return apt, bpt, ant, bnt
end
