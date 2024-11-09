# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Spheroidal Wave Functions

pro_rad1(m, n, c, x),
pro_rad2(m, n, c, x)
+ specfun_segv, specfun_rswfp

pro_ang1(m, n, c, x),
obl_ang1(m, n, c, x)
+ specfun_segv, specfun_aswfa

obl_rad1(m, n, c, x),
obl_rad2(m, n, c, x)
+ specfun_segv, specfun_rswfo

pro_cv(m, n, c),
obl_cv(m, n, c),
pro_cv_seq(m, n, c),
obl_cv_seq(m, n, c)
+ specfun_segv
"""

"""

Compute the characteristic values of spheroidal wave functions

Input:  
m  --- Mode parameter
n  --- Mode parameter
c  --- Spheroidal parameter
KD --- Function code
    KD=1 for Prolate; KD=-1 for Oblate

Output:  
CV --- Characteristic value for given m, n and c
EG(L) --- Characteristic value for mode m and n'
        ( L = n' - m + 1 )
"""
function segv(m::Int, n::Int, c::T, kd::Int, eg::Vector{T}) where {T<:AbstractFloat}
    @assert kd==1 || kd==-1 "Bad kd, kd not in [1, -1]"
    eg_len = n - m + 1
    @assert eg_len > 0 "Bad (m, n)"
    @assert length(eg) >= 200 "eg[] too small, need length(eg) >= 200"

    if c < T(1e-10)
        for i in 1:eg_len
            eg[i] = (i + m) * (i + m - 1)
        end
        cv = eg[eg_len]
        return cv, eg
    end

    # TODO-opt: Following array sizes should be decided dynamically
    a = fill(T(0), 300)
    b = fill(T(0), 100)
    cv0 = fill(T(0), 100)
    d = fill(T(0), 300)
    e = fill(T(0), 300)
    f = fill(T(0), 300)
    g = fill(T(0), 300)
    h = fill(T(0), 100)

    icm = Int((n - m + 2) ÷ 2)
    nm = Int(10 + trunc(Int, 0.5 * (n - m) + c))
    cs = T(c * c * kd)
    k = Int(0)
    @assert nm <= 300 "a[], d[], g[] out of range"
    @assert nm <= 300 "e[], k[] out of range"
    @assert icm <= 100 "h[] out of range"
    for l in 0:1
        for i in 1:nm
            k = (l == 0 ? 2 * (i-1) : 2*i-1)
            dk0 = m + k
            dk1 = m + k + 1
            dk2 = 2 * (m + k)
            d2k = 2 * m + k
            a[i] = ((d2k+2)*(d2k+1)) / ((dk2+3)*(dk2+5)) * cs
            d[i] = dk0*dk1 + (2*dk0*dk1 - 2*m*m - 1) / ((dk2-1)*(dk2+3)) * cs
            g[i] = k*(k-1) / ((dk2-3)*(dk2-1)) * cs
        end
        for k in 2:nm
            e[k] = sqrt(a[k-1] * g[k])
            f[k] = e[k] * e[k]
        end
        f[1] = T(0)
        e[1] = T(0)

        xa = d[nm] + abs(e[nm])
        xb = d[nm] - abs(e[nm])
        nm1 = Int(nm - 1)
        @assert nm1 <= 300 "d[] out of range"
        for i in 1:nm1
            t = abs(e[i]) + abs(e[i+1])
            t1 = d[i] + t
            xa = max(xa, t1)
            t1 = d[i] - t
            xb = min(xb, t1)
        end
        for i in 1:icm
            b[i] = xa
            h[i] = xb
        end

        @assert (2*icm) <= 200 "eg[2*k] out of range"
        for k in 1:icm
            for k1 in k:icm
                if b[k1] < b[k]
                    b[k] = b[k1]
                    break
                end
            end

            if k != 1 && h[k] < h[k-1]
                h[k] = h[k-1]
            end

            x1 = T(0.0)
            while true
                x1 = (b[k] + h[k]) / T(2.0)
                cv0[k] = x1
                if abs((b[k] - h[k]) / x1) < 1e-14
                    break
                end

                j = Int(0)
                s = T(1.0)
                for i in 1:nm
                    if s == 0.0
                        s += T(1e-30)
                    end
                    t = f[i] / s
                    s = d[i] - t - x1
                    if s < 0.0
                        j += 1
                    end
                end
                if j < k
                    h[k] = x1
                else
                    b[k] = x1
                    if j >= icm
                        b[icm] = x1
                    else
                        h[j+1] = max(h[j+1], x1)
                        b[j] = min(b[j], x1)
                    end
                end
            end

            cv0[k] = x1
            if l == 0
                eg[2*k - 1] = cv0[k]
            elseif l == 1
                eg[2*k] = cv0[k]
            end
        end
    end

    cv = eg[eg_len]
    return cv, eg
end
