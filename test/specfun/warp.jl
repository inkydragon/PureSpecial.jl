# SPDX-License-Identifier: MIT
using Libdl

"""
## Compile fortran lib
```sh
make -C test/f77/
```
"""
# where to find: `specfun.so`
github_workspace = get(ENV, "GITHUB_WORKSPACE", nothing)
if isnothing(github_workspace)
    @warn """Cannot find f77 ref impl lib!  You may want to set `ENV["GITHUB_WORKSPACE"]`"""
else
    github_workspace = joinpath(github_workspace, "test/f77/")
    @info "Load ref impl for: '$(github_workspace)'"
    push!(DL_LOAD_PATH, github_workspace)
end
const libspecfun = Libdl.dlopen("libspecfun.so")


"Loading symbols from libspecfun."
f77func(s::Symbol) = Libdl.dlsym(libspecfun, Symbol("$(s)_")) 


#= ## Airy functions =#
"""
Warp fortran `specfun.AIRYB`.

## Output
- `(ai, bi, ad, bd) :: NTuple{4, Float64}`
"""
function _airyb(z::Float64)
    ai, bi, ad, bd = Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN)
    # AIRYB(X,AI,BI,AD,BD)
    # void specfun_airyb(double, double *, double *, double *, double *);
    ret = ccall(f77func(:airyb), Cvoid,
        (Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        z,
        ai, bi, ad, bd)

    ai[], bi[], ad[], bd[]
end

"""
Warp fortran `specfun.AIRYZO`.

- Input: `nt`, `kf`
- Output: `xa, xb, xc, xd`
"""
function _airyzo!(
    nt::Int, kf::Int, 
    xa::Vector{Float64}, xb::Vector{Float64},
    xc::Vector{Float64}, xd::Vector{Float64})
    # AIRYZO(NT,KF,XA,XB,XC,XD)
    # void specfun_airyzo(int nt, int kf, double *xa, double *xb, double *xc, double *xd);
    ret = ccall(f77func(:airyzo), Cvoid,
        (Ref{Int64}, Ref{Int64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        nt, kf,
        xa, xb, xc, xd)
end

"""
Warp fortran `specfun.ITAIRY`.

- Input: `x`
- Output: `apt, bpt, ant, bnt`
"""
function _itairy(x::Float64)
    @assert x >= 0
    apt, bpt, ant, bnt =
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN)
    # SUBROUTINE ITAIRY(X,APT,BPT,ANT,BNT)
    # void specfun_itairy(double x,
    #   double *apt, double *bpt, double *ant, double *bnt);
    ccall(f77func(:itairy), Cvoid,
        (Ref{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        x,
        apt, bpt, ant, bnt)
    apt[], bpt[], ant[], bnt[]
end


#= ## Bessel functions =#
"""
Warp fortran `specfun.GAM0`.
- Input: `x`
- Output: `ga`
"""
function _gam0(x::Float64)
    ga = Ref{Float64}(NaN)
    # SUBROUTINE GAM0 (X,GA)
    # double specfun_gam0(double x);
    ccall(f77func(:gam0), Cvoid,
        (Ref{Float64}, Ref{Float64}),
        x, ga)

    ga[]
end

"""
Warp fortran `specfun.GAMMA2`.
- Input: `x`
- Output: `ga`
"""
function _gamma2(x::Float64)
    ga = Ref{Float64}(NaN)
    # SUBROUTINE GAMMA2(X,GA)
    # double specfun_gamma2(double x);
    ccall(f77func(:gamma2), Cvoid,
        (Ref{Float64}, Ref{Float64}),
        x, ga)

    ga[]
end


#= ## Zeros of Bessel functions =#
"""
Warp fortran `specfun.BJNDD`.

- Input: `x`, `n`
- Output: `bj[], dj[], fj[]`

!!! note
    Fortran code uses fixed-size array.
"""
function _bjndd(x::Float64, n::Int)
    # Pass fixed-size array.
    bj, dj, fj = zeros(Float64, 101), zeros(Float64, 101), zeros(Float64, 101)
    # BJNDD(N,X,BJ,DJ,FJ)
    # void specfun_bjndd(double x, int n, double *bj, double *dj, double *fj);
    ret = ccall(f77func(:bjndd), Cvoid,
        (Ref{Int64}, Ref{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        n, x,
        bj, dj, fj)

    bj, dj, fj
end

"""
Warp fortran `specfun.JDZO`.

- Input: `nt`
- Output: `zo, n, m, p`
"""
function _jdzo(nt::Int)
    @assert nt in 1:1200
    zo = zeros(Float64, 1401)
    n = zeros(Int32, 1400)
    m = zeros(Int32, 1400)
    p = zeros(Int32, 1400)
    # NOTE: different params orders.
    #   JDZO(NT,N,M,P,ZO)
    #   void specfun_jdzo(int nt, double *zo, int *n, int *m, int *p);
    ccall(f77func(:jdzo), Cvoid,
        (Ref{Int32},
         Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        nt,
        n, m, p, zo)

    zo[1:nt], n[1:nt], m[1:nt], p[1:nt]
end

"""
Warp fortran `specfun.MSTA1`.
"""
function _msta1(x::Float64, mp::Int32)
    # INTEGER FUNCTION MSTA1(X,MP)
    # int specfun_msta1(double x, int mp);
    ccall(f77func(:msta1), Cint,
        (Ref{Float64}, Ref{Int32}),
        x, mp)
end
_msta1(x::Float64, mp::Int64) = _msta1(x, Int32(mp))

"""
Warp fortran `specfun.MSTA2`.
"""
function _msta2(x::Float64, n::Int32, mp::Int32)
    # INTEGER FUNCTION MSTA2(X,N,MP)
    # int specfun_msta2(double x, int n, int mp);
    ccall(f77func(:msta2), Cint,
        (Ref{Float64}, Ref{Int32}, Ref{Int32}),
        x, n, mp)
end
_msta2(x::Float64, n::Int64, mp::Int64) = _msta2(x, Int32(n), Int32(mp))

"""
Warp fortran `specfun.JYNBH`.

- Output: `bj ,by`
- Return: `nm`
"""
function _jynbh!(
    x::Float64, n::Int32, nmin::Int32,
    bj::Vector{Float64}, by::Vector{Float64})
    nm = Ref{Int32}(0)
    # SUBROUTINE JYNBH(N,NMIN,X, NM,BJ,BY)
    # void specfun_jynbh(int n, int nmin, double x, 
    #   int *nm, double *bj, double *by);
    ccall(f77func(:jynbh), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64},
         Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        n, nmin, x,
        nm, bj, by)

    nm[]
end
_jynbh!(x, n::Int64, nmin::Int64, bj, by) =
    _jynbh!(x, Int32(n), Int32(nmin), bj, by)

"""
Warp fortran `specfun.JYNDD`.

- Return: `bjn, djn, fjn, byn, dyn, fyn`
"""
function _jyndd(x::Float64, n::Int32)
    @assert x > 0.0

    bjn, djn, fjn, byn, dyn, fyn =
        Ref{Float64}(0), Ref{Float64}(0), Ref{Float64}(0),
        Ref{Float64}(0), Ref{Float64}(0), Ref{Float64}(0)
    # SUBROUTINE JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
    # void specfun_jyndd(int n, double x,
    #   double *bjn, double *djn, double *fjn,
    #   double *byn, double *dyn, double *fyn);
    ccall(f77func(:jyndd), Cvoid,
        (Ref{Int32}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}),
        n, x,
        bjn, djn, fjn,
        byn, dyn, fyn)

    bjn[], djn[], fjn[], byn[], dyn[], fyn[]
end
_jyndd(x, n::Int64) = _jyndd(x, Int32(n))

"""
Warp fortran `specfun.JYZO`.

- Input: `n, nt`
- Output: `rj0, rj1, ry0, ry1`
"""
function _jyzo!(n::Int64, nt::Int64,
    rj0::Vector{Float64}, rj1::Vector{Float64},
    ry0::Vector{Float64}, ry1::Vector{Float64})
    @assert n >= 0
    @assert length(rj0) >= nt
    @assert length(rj1) >= nt
    @assert length(ry0) >= nt
    @assert length(ry1) >= nt
    #   SUBROUTINE JYZO(N,NT,RJ0,RJ1,RY0,RY1)
    #   void specfun_jyzo(int n, int nt,
    #       double *rj0, double *rj1, double *ry0, double *ry1);
    ccall(f77func(:jyzo), Cvoid,
        (Ref{Int32}, Ref{Int32},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        n, nt,
        rj0, rj1, ry0, ry1)
end


"""Parabolic Cylinder functions

- [ ] _pbwa
    - [x] _cgama

- [ ] _pbvv
    - [ ] _gamma2
    - [ ] _vvla
    - [ ] _vvsa

- [ ] _pbdv
    - [ ] _dvsa
    - [ ] _dvla

"""


#= ## Kelvin functions =#
"""
Warp fortran `specfun.KLVNA`.

- Input: `x`
- Output: `ber, bei, ger, gei, der, dei, her, hei`
"""
function _klvna(x::Float64)
    ber, bei, ger, gei, der, dei, her, hei = 
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN),
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN) 
    # KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
    # void specfun_klvna(double x, 
    #     double *ber, double *bei, double *ger, double *gei,
    #     double *der, double *dei, double *her, double *hei);
    ret = ccall(f77func(:klvna), Cvoid,
        (Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        x,
        ber, bei, ger, gei,
        der, dei, her, hei)
        
    ber[], bei[], ger[], gei[], der[], dei[], her[], hei[]
end

"""
Warp fortran `specfun.KLVNZO`.

- Input: `nt`, `kd`
- Output: `zo::Vector{Float64}`
"""
function _klvnzo(nt::Int64, kd::Int64)
    zo = zeros(Float64, nt)
    # KLVNZO(NT,KD,ZO)
    # void specfun_klvnzo(int nt, int kd, double *zo);
    ret = ccall(f77func(:klvnzo), Cvoid,
        (Ref{Int64}, Ref{Int64},
         Ptr{Float64}),
        nt, kd,
        zo)

    zo
end


"""
Warp fortran `specfun.CGAMA`.

- Input: `X,Y, KF`
- Output: `GR,GI`
"""
function _cgama(z::Complex{Float64}, kf::Int)
    gr = Ref{Float64}(NaN)
    gi = Ref{Float64}(NaN)
    # SUBROUTINE CGAMA(X,Y,KF,GR,GI)
    # double complex specfun_cgama(double complex z, int kf);
    ccall(f77func(:cgama), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Int32},
         Ref{Float64}, Ref{Float64}),
        real(z), imag(z), Int32(kf),
        gr, gi)

    complex(gr[], gi[])
end

"""
Warp fortran `specfun.CCHG`.

- Input: `A,B,Z`
- Output: `CHG`
"""
function _cchg(a::Float64, b::Float64, z::Complex{Float64})
    chg = Ref{Complex{Float64}}(NaN + NaN*im)
    # SUBROUTINE CCHG(A,B,Z,CHG)
    # double complex specfun_cchg(double a, double b, double complex z);
    ccall(f77func(:cchg), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Complex{Float64}},
         Ref{Complex{Float64}}),
        a, b, z,
        chg)

    chg[]
end

"""
Warp fortran `specfun.CHGM`.

- Input: `A,B,X`
- Output: `HG`
"""
function _chgm(a::Float64, b::Float64, x::Float64)
    hg = Ref{Float64}(NaN)
    # SUBROUTINE CHGM(A,B,X,HG)
    # double specfun_chgm(double x, double a, double b);
    ccall(f77func(:chgm), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}),
        a, b, x,
        hg)

    hg[]
end
