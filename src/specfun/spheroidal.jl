# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Spheroidal Wave Functions

pro_rad1(m, n, c, x),
pro_rad2(m, n, c, x)
+ specfun_segv, specfun_rswfp
    -> sdmn
    -> rmn1 -> sckb,sphj
        -> msta1,msta2
    -> rmn2l -> sphy
    -> rmn2sp -> kmn,lqmns,lpmns

pro_ang1(m, n, c, x),
obl_ang1(m, n, c, x)
+ specfun_segv, specfun_aswfa
    -> sdmn,sckb

obl_rad1(m, n, c, x),
obl_rad2(m, n, c, x)
+ specfun_segv, specfun_rswfo
    -> sdmn,sckb,[rmn2l],rmn2so
        -> qstar,cbk,gmn,kmn,[rmn1]

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

"""
Compute the expansion coefficients of the
prolate and oblate spheroidal functions, dk

Input: 
m  --- Mode parameter
n  --- Mode parameter
c  --- Spheroidal parameter
cv --- Characteristic value
KD --- Function code
    KD=1 for prolate; KD=-1 for oblate

Output: 
DF(k) --- Expansion coefficients dk;
    DF(1), DF(2), ... correspond to
    d0, d2, ... for even n-m and d1,
    d3, ... for odd n-m
"""
function sdmn!(m::Int, n::Int, c::T, cv::T, kd::Int, df::Vector{T}) where {T<:AbstractFloat}
    @assert kd in Set([1, -1]) "Bad kd: kd not in [1, -1]"
    @assert length(df) >= 200 "df[] too small, need length(df) >= 200"

    nm = 25 + trunc(Int, 0.5*(n-m) + c)
    @assert nm <= length(df) "(n-m) too large"
    @assert ((n-m)÷2 + 1) <= length(df) "((n-m)÷2 + 1) too large"
    if c < 1e-10
        for i in 1:nm
            df[i] = 0
        end
        df[(n-m)÷2 + 1] = 1
        return
    end

    a = fill(T(0), nm + 2)
    d = fill(T(0), nm + 2)
    g = fill(T(0), nm + 2)
    cs = c*c*kd
    ip = (n - m) % 2

    for i in 1:(nm+2)
        k = (ip == 0 ? 2*(i-1) : 2*i - 1)
        dk0 = m + k
        dk1 = m + k + 1
        dk2 = 2 * (m + k)
        d2k = 2 * m + k
        a[i] = (d2k+2)*(d2k+1) / ((dk2+3)*(dk2+5)) * cs
        d[i] = dk0*dk1 + (2*dk0*dk1 - 2*m*m - 1) / ((dk2-1)*(dk2+3)) * cs
        g[i] = k * (k-1) / ((dk2-3)*(dk2-1)) * cs
    end

    fs = T(1.0)
    f1 = T(0.0)
    f0 = T(1e-100)
    kb = Int(0)
    @assert (nm+1) <= length(df) "(n-m) too large"
    df[nm+1] = T(0.0)
    fl = T(0.0)
    for k in nm:-1:1
        f = -((d[k+1]-cv)*f0 + a[k+1]*f1) / g[k+1]
        if abs(f) > abs(df[k+1])
            df[k] = f
            f1 = f0
            f0 = f

            if abs(f) > T(1e100)
                for k1 in k:nm
                    df[k1] *= T(1e-100)
                end
                f1 *= T(1e-100)
                f0 *= T(1e-100)
            end
        else
            kb = k
            fl = df[k+1]
            f1 = T(1e-100)
            f2 = -((d[1]-cv) / a[1]) * f1
            df[1] = f1

            if kb == 1
                fs = f2
            elseif kb == 2
                df[2] = f2
                fs = -((d[2]-cv)*f2 + g[2]*f1) / a[2]
            else
                df[2] = f2
                for j in 3:(kb+1)
                    f = -((d[j-1]-cv)*f2 + g[j-1]*f1) / a[j-1]
                    if j <= kb
                        df[j] = f
                    end
                    if abs(f) > T(1e100)
                        for k1 in 1:j
                            df[k1] *= T(1e-100)
                        end
                        f *= T(1e-100)
                        f2 *= T(1e-100)
                    end
                    f1 = f2
                    f2 = f
                end
                fs = f
            end

            break
        end
    end

    su1 = T(0.0)
    r1 = T(1.0)
    for j in range(m+ip+1, 2*(m+ip))
        r1 *= j
    end

    su1 = df[1] * r1
    for k in 2:kb
        r1 = -r1 * (k+m+ip-1.5) / (k-1.0)
        su1 += r1 * df[k]
    end

    su2 = T(0.0)
    sw = T(0.0)
    for k in (kb+1):nm
        if k != 1
            r1 = -r1 * (k+m+ip-1.5) / (k-1.0)
        end
        su2 += r1 * df[k]
        if abs(sw - su2) < abs(su2)*T(1e-14)
            break
        end

        sw = su2
    end

    r3 = T(1.0)
    for j in 1:((m+n+ip) ÷ 2)
        r3 *= j + 0.5*(m+n+ip)
    end
    r4 = T(1.0)
    for j in 1:((n-m-ip) ÷ 2)
        r4 *= -4.0 * j
    end

    s0 = r3 / (fl*(su1/fs) + su2) / r4
    for k in 1:kb
        df[k] *= fl / fs * s0
    end
    for k in (kb+1):nm
        df[k] *= s0
    end
end

"""
Compute the expansion coefficients of the
prolate and oblate spheroidal functions

Input:
m  --- Mode parameter
n  --- Mode parameter
c  --- Spheroidal parameter
DF(k) --- Expansion coefficients dk

Output:  
CK(k) --- Expansion coefficients ck;
    CK(1), CK(2), ... correspond to
    c0, c2, ...
"""
function sckb!(m::Int, n::Int, c::T, df::Vector{T}, ck::Vector{T}) where {T<:AbstractFloat}
    if c <= T(1.0e-10)
        c = T(1.0e-10)
    end
    nm = 25 + trunc(Int, 0.5*(n-m) + c)
    ip = (n - m) % 2
    reg = ((m+nm) > 80 ? T(1.0e-200) : T(1.0))
    fac = -T(0.5)^m
    sw = T(0.0)
    for k in 0:(nm-1)
        fac = -fac
        i1 = 2*k + ip + 1
        r = reg
        for i in i1:(i1 + 2*m - 1)
            r *= i
        end
        i2 = k + m + ip
        for i in i2:(i2 + k - 1)
            r *= i + T(0.5)
        end

        sum = r * df[k+1]
        for i in (k+1):nm
            d1 = T(2.0)*i + ip
            d2 = T(2.0)*m + d1
            d3 = i + m + ip - T(0.5)
            r = r*d2*(d2-1)*i*(d3+k) / (d1*(d1-1)*(i-k)*d3)
            sum += r * df[i+1]
            if abs(sw - sum) < abs(sum) * T(1.0e-14)
                break
            end

            sw = sum
        end
        r1 = reg
        for i in 2:(m+k)
            r1 *= i
        end
        ck[k+1] = fac * sum / r1
    end
end

"""
Compute the prolate and oblate spheroidal angular
functions of the first kind and their derivatives

Input :  
m  --- Mode parameter,  m = 0,1,2,...
n  --- Mode parameter,  n = m,m+1,...
c  --- Spheroidal parameter
x  --- Argument of angular function, |x| ≤ 1.0
KD --- Function code
    KD=1 for prolate;  KD=-1 for oblate
cv --- Characteristic value

Output:  
s1f --- Angular function of the first kind
s1d --- Derivative of the angular function of
        the first kind

Call: sdmn,sckb
"""
function aswfa(m::Int, n::Int, c::T, x::T, kd::Int, cv::T) where {T<:AbstractFloat}
    @assert abs(x) <= 1 "Bad x: |x| ≤ 1.0"
    _EPS = T(1e-14)
    ip = ((n - m) % 2 == 0) ? 0 : 1
    nm = trunc(Int, 40 + (n - m) / 2 + c)
    nm2 = nm ÷ 2 - 2
    x0 = x
    x = abs(x)

    df = zeros(T, 200)
    ck = zeros(T, 200)
    sdmn!(m, n, c, cv, kd, df)
    sckb!(m, n, c, df, ck)

    x1 = T(1) - x * x
    a0 = ((m == 0) && (x1 == 0.0)) ? T(1.0) : x1^(T(0.5) * m)

    su1 = ck[1]
    for k in 1:nm2
        r = ck[k+1] * x1^k
        su1 += r
        if (k >= 10) && (abs(r/su1) < _EPS)
            break
        end
    end

    s1f = a0 * x^ip * su1
    s1d = T(0.0)
    if x == T(1.0)
        if m == 0
            s1d = ip*ck[1] - 2.0*ck[2]
        elseif m == 1
            s1d = -1e100
        elseif m == 2
            s1d = -2.0*ck[1]
        elseif m == 3
            s1d = 0.0
        end
    else
        d0 = ip - m / x1 * x^(ip+1)
        d1 = -2.0 * a0 * x^(ip+1)
        su2 = ck[2]
        for k in 2:nm2
            r = k * ck[k+1] * x1^(k-1)
            su2 += r
            if (k >= 10) && (abs(r/su2) < _EPS)
                break
            end
        end
        s1d = d0 * a0 * su1 + d1 * su2
    end

    if (x0 < 0.0) && (ip == 0)
        s1d = -s1d
    end
    if (x0 < 0.0) && (ip == 1)
        s1f = -s1f
    end

    return s1f, s1d
end


"""

Compute spherical Bessel functions jn(x) and their derivatives

Input:  
x --- Argument of jn(x)
n --- Order of jn(x)  ( n = 0,1,… )

Output:  
sj --- Array of jn(x) values
dj --- Array of jn'(x) values
nm --- Highest order computed

Routines called:
- msta1
- msta2
"""
function sphj!(n::Int, x::T, sj::Vector{Float64}, dj::Vector{Float64}) where {T<:AbstractFloat}
    @assert n >= 0
    # NOTE: in f77 version
    #   DIMENSION SJ(0:N),DJ(0:N)
    @assert length(sj) >= 1
    @assert length(dj) >= 1

    nm = n
    if abs(x) < T(1e-100)
        for k in 0:n
            sj[k+1] = T(0.0)
            dj[k+1] = T(0.0)
        end
        sj[0+1] = T(1.0)
        if n > 0
            dj[1+1] = T(1.0) / 3.0
        end
        return sj, dj, nm
    end

    sj[0+1] = sin(x) / x
    dj[0+1] = (cos(x) - sin(x) / x) / x
    if n < 1
        return sj, dj, nm
    end

    sj[1+1] = (sj[0+1] - cos(x)) / x
    if n >= 2
        sa = sj[0+1]
        sb = sj[1+1]
        m = msta1(x, 200)
        if m < n
            nm = m
        else
            m = msta2(x, n, 15)
        end

        f = T(0.0)
        f0 = T(0.0)
        f1 = T(1e-100)
        for k in m:-1:0
            f = ((2*k+3) * f1 / x) - f0
            if k <= nm
                sj[k+1] = f
            end
            f0 = f1
            f1 = f
        end

        cs = abs(sa) > abs(sb) ? sa/f : sb/f0
        for k in 0:nm
            sj[k+1] *= cs
        end
    end
    for k in 1:nm
        dj[k+1] = sj[k] - (k+1)*sj[k+1] / x
    end

    return sj, dj, nm
end

"""
Compute prolate and oblate spheroidal radial
functions of the first kind for given m, n, c, and x

Routines called:
- sckb
- sphj
"""
function rmn1(m::Int, n::Int, c::T, x::T, kd::Int, df::Vector{T}) where {T<:AbstractFloat}
    @assert (1.0 - kd / (x * x)) >= 0 
    _EPS = 1e-14
    ck = zeros(T, 200)

    nm1 = trunc(Int, (n - m) / 2)
    ip = ifelse((n-m) == 2*nm1, 0, 1)
    nm = 25 + nm1 + trunc(Int, c)
    reg = ifelse((m+nm) > 80, T(1e-200), T(1.0))
    r0 = reg
    for j in 1:(2*m + ip)
        r0 *= j
    end

    r = T(r0)
    suc = r * df[1]
    sw = T(0.0)
    for k in 2:nm
        r *= (m + k - 1) * (m + k + ip - 1.5) / ((k - 1) * (k + ip - 1.5))
        suc += r * df[k]
        if (k > nm1) && (abs(suc - sw) < abs(suc) * _EPS)
            break
        end
    
        sw = suc
    end

    if x == 0.0
        sckb!(m, n, c, df, ck)

        sum = T(0.0)
        sw1 = T(0.0)
        for j in 1:nm
            sum += ck[j]
            if abs(sum - sw1) < abs(sum) * _EPS
                break
            end

            sw1 = sum
        end

        r1 = T(1.0)
        for j in 1:trunc(Int, (n + m + ip) / 2)
            r1 *= j + 0.5 * (n + m + ip)
        end
        r2 = T(1.0)
        for j in 1:m
            r2 *= 2.0 * c * j
        end
        r3 = T(1.0)
        for j in 1:trunc(Int, (n - m - ip) / 2)
            r3 *= j
        end
        sa0 = (2*(m+ip) + 1) * r1 / (2.0^n * c^ip * r2 * r3)
        if ip == 0
            r1f = sum / (sa0 * suc) * df[1] * reg
            r1d = T(0.0)
        else
            r1f = T(0.0)
            r1d = sum / (sa0 * suc) * df[1] * reg
        end

        return r1f, r1d
    end

    # NOTE: in f77 version
    #   DIMENSION SJ(0:251),DJ(0:251)
    sj = zeros(T, 251+1)
    dj = zeros(T, 251+1)
    @assert ((m + 2*nm - 2 + ip) + 1) <= (251+1) "sj[], dj[] Out of index"

    cx = c * x
    nm2 = 2 * nm + m
    _, _, nm2 = sphj!(nm2, cx, sj, dj)

    a0 = (1.0 - kd / (x * x))^(0.5 * m) / suc
    r1f = T(0.0)
    sw = T(0.0)
    for k in 1:nm
        l = 2*k + m - n - 2 + ip
        lg = ifelse(l % 4 == 0, 1, -1)
        if k == 1
            r = r0
        else
            r *= (m + k - 1.0) * (m + k + ip - 1.5) / ((k - 1.0) * (k + ip - 1.5))
        end
        np = m + 2 * k - 2 + ip
        r1f += lg * r * df[k] * sj[np + 1]
        if (k > nm1) && (abs(r1f - sw) < abs(r1f) * _EPS)
            break
        end

        sw = r1f
    end

    r1f *= a0
    b0 = kd * m / (x^3 * (1.0 - kd / (x * x))) * r1f
    sud = T(0.0)
    sw = T(0.0)
    for k in 1:nm
        l = 2 * k + m - n - 2 + ip
        lg = ifelse(l % 4 == 0, 1, -1)
        if k == 1
            r = r0
        else
            r *= (m + k - 1) * (m + k + ip - 1.5) / ((k-1) * (k + ip - 1.5))
        end
        np = m + 2 * k - 2 + ip
        sud += lg * r * df[k] * dj[np + 1]
        if (k > nm1) && (abs(sud - sw) < abs(sud) * _EPS)
            break
        end

        sw = sud
    end

    r1d = b0 + a0 * c * sud
    return r1f, r1d
end

"""
Compute spherical Bessel functions yn(x) and
their derivatives

Input :  
x --- Argument of yn(x) ( x ≥ 0 )
n --- Order of yn(x) ( n = 0,1,… )

Output:  
SY(n) --- yn(x)
DY(n) --- yn'(x)
NM --- Highest order computed
"""
function sphy!(n::Int, x::T, sy::Vector{Float64}, dy::Vector{Float64}) where {T<:AbstractFloat}
    @assert n >= 0
    @assert x >= 0
    @assert length(sy) >= (n+1)
    @assert length(dy) >= (n+1)

    nm = n
    if x < 1.0e-60
        fill_len = n + 1
        sy_view = @view sy[1:fill_len]
        fill!(sy_view, -1.0e300)
        dy_view = @view dy[1:fill_len]
        fill!(dy_view, 1.0e300)
        return sy, dy, nm
    end

    sy[0+1] = -cos(x) / x
    f0 = sy[0+1]
    dy[0+1] = (sin(x) + cos(x) / x) / x
    if n < 1
        return sy, dy, nm
    end

    sy[1+1] = (sy[0+1] - sin(x)) / x
    f1 = sy[1+1]
    for k in 2:n
        f = (2.0*k - 1.0) * f1 / x - f0
        sy[k+1] = f
        if abs(f) >= 1.0e300
            nm = k - 1
            return sy, dy, nm
        end

        f0 = f1
        f1 = f
    end
    for k in 1:nm
        dy[k+1] = sy[k] - (k + 1.0) * sy[k+1] / x
    end

    return sy, dy, nm
end

"""
Compute prolate and oblate spheroidal radial
functions of the second kind for given m, n,
c and a large cx

Routine called:
SPHY for computing the spherical Bessel
functions of the second kind
"""
function rmn2l(m::Int, n::Int, c::T, x::T, kd::Int, df::Vector{T}) where {T<:AbstractFloat}
    @assert length(df) >= 200
    _EPS = 1.0e-14
    # In f77 version:
    #   DIMENSION DF(200),SY(0:251),DY(0:251)
    sy = zeros(T, 251+1)
    dy = zeros(T, 251+1)

    ip = ifelse((n - m) % 2 == 0, 0, 1)
    nm1 = (n - m) ÷ 2
    nm = 25 + nm1 + trunc(Int, c)
    reg = ifelse(m + nm > 80, 1.0e-200, 1.0)
    nm2 = 2 * nm + m
    cx = c * x
    _, _, nm2 = sphy!(nm2, cx, sy, dy)
    r0 = reg
    for j in 1:(2*m + ip)
        r0 *= j
    end

    r = r0
    suc = r * df[1]
    sw = 0.0
    for k in 2:nm
        r *= (m + k - 1.0) * (m + k + ip - 1.5) / ((k - 1.0) * (k + ip - 1.5))
        suc += r * df[k]
        if k > nm1 && abs(suc - sw) < abs(suc) * _EPS
            break
        end

        sw = suc
    end

    a0 = (1.0 - kd / (x * x))^(0.5 * m) / suc
    r2f = 0.0
    r2d = 0.0
    eps1 = 0.0
    np = 0
    for k in 1:nm
        l = 2 * k + m - n - 2 + ip
        lg = ifelse(l % 4 == 0, 1, -1)
        if k == 1
            r = r0
        else
            r *= (m + k - 1.0) * (m + k + ip - 1.5) / ((k - 1.0) * (k + ip - 1.5))
        end
        np = m + 2 * k - 2 + ip
        r2f += lg * r * (df[k] * sy[np + 1])
        eps1 = abs(r2f - sw)
        if k > nm1 && eps1 < abs(r2f) * _EPS
            break
        end

        sw = r2f
    end

    id1 = trunc(Int, log10(eps1 / abs(r2f) + _EPS))
    r2f *= a0
    if np >= nm2
        id = 10
        return T(r2f), T(r2d), id
    end

    b0 = kd * m / x^3.0 / (1.0 - kd / (x * x)) * r2f
    sud = 0.0
    eps2 = 0.0
    for k in 1:nm
        l = 2 * k + m - n - 2 + ip
        lg = ifelse(l % 4 == 0, 1, -1)
        if k == 1
            r = r0
        else
            r *= (m + k - 1.0) * (m + k + ip - 1.5) / ((k - 1.0) * (k + ip - 1.5))
        end
        np = m + 2 * k - 2 + ip
        sud += lg * r * (df[k] * dy[np + 1])
        eps2 = abs(sud - sw)
        if k > nm1 && eps2 < abs(sud) * _EPS
            break
        end

        sw = sud
    end

    r2d = b0 + a0 * c * sud
    id2 = trunc(Int, log10(eps2 / abs(sud) + _EPS))
    id = max(id1, id2)
    return T(r2f), T(r2d), id
end

"""
Compute prolate spheriodal radial functions of the
first and second kinds, and their derivatives

Input:  
m  --- Mode parameter, m = 0,1,2,...
n  --- Mode parameter, n = m,m+1,m+2,...
c  --- Spheroidal parameter
x  --- Argument of radial function ( x > 1.0 )
cv --- Characteristic value
KF --- Function code
    KF=1 for the first kind
    KF=2 for the second kind
    KF=3 for both the first and second kinds

Output:  
R1F --- Radial function of the first kind
R1D --- Derivative of the radial function of
        the first kind
R2F --- Radial function of the second kind
R2D --- Derivative of the radial function of
        the second kind

Routines called:
- sdmn
- rmn1
- rmn2l
- rmn2sp
"""
function rswfp(m::Int, n::Int, c::T, x::T, cv::T, kf::Int) where {T<:AbstractFloat}
    df = zeros(T, 200)
    kd = 1

    sdmn!(m, n, c, cv, kd, df)

    r1f = 0.0
    r1d = 0.0
    r2f = 0.0
    r2d = 0.0
    if kf != 2
        r1f, r1d = rmn1(m, n, c, x, kd, df)
    end
    # if kf > 1
    #     r2f, r2d, id = rmn2l(m, n, c, x, kd, df)
    #     if id > -8
    #         r2f, r2d = rmn2sp(m, n, c, x, cv, kd, df)
    #     end
    # end

    # return r1f, r1d, r2f, r2d
end
