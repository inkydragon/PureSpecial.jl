# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md 
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
    # TODO: replace by pi
    pi64 = 3.141592653589793
    @assert Float64(pi) === pi64
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
            if abs(r) < (abs(fx) * eps)
                break
            end
        end #= 10 =#

        #= 15 =#
        gx = x
        r = x
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k + 1.0) * x
            gx += r
            if abs(r) < (abs(gx) * eps)
                break
            end
        end #= 20 =#
        #= 25 =#
        ai = c1 * fx - c2 * gx
        bi = sr3 * (c1 * fx + c2 * gx)

        df = 0.5 * x * x
        r = df
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k + 2.0) * x
            df += r
            if abs(r) < (abs(df) * eps)
                break
            end
        end #= 30 =#

        #= 35 =#
        dg = 1.0
        r = 1.0
        for k = 1:40
            r *= x / (3.0 * k) * x / (3.0 * k - 2.0) * x
            dg += r
            if abs(r) < (abs(dg) * eps)
                break
            end
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

            xcs = cos(xe + pi64 / 4.0)
            xss = sin(xe + pi64 / 4.0)

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
