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


#= gamma.jl 
    + gam0
    + gamma2
    + gaih
    + cgama
    + psi (PSI_SPEC)
=#
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

"""
Warp fortran `specfun.GAIH`.
- Input: `x`
- Output: `ga`
"""
function _gaih(x::Float64)
    ga = Ref{Float64}(NaN)
    # SUBROUTINE GAIH(X,GA)
    # double specfun_gaih(double x);
    ccall(f77func(:gaih), Cvoid,
        (Ref{Float64}, Ref{Float64}),
        x, ga)

    ga[]
end

"""
Warp fortran `specfun.CGAMA`.
- Input: `z, kf`
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
    SUBROUTINE PSI_SPEC(X,  PS)
    double psi_spec(double x)

- Input: `x`
- Output: `psi(x)`
"""
function _psi(x::Float64)
    psi = Ref{Float64}(NaN)
    ccall(f77func(:psi_spec), Cvoid,
        (Ref{Float64}, Ref{Float64}),
        x, psi)
    psi[]
end


"""Exponential and Trigonometric Integrals"""
#=
- e1xb
- [ ] e1xa
- [ ] enxb
- [ ] enxa
- e1z
- eix
- eixz
=#

function _e1xb(x::Float64)
    e1 = Ref{Float64}(NaN)
    # SUBROUTINE E1XB(X,E1)
    # double specfun_e1xb(double x);
    ccall(f77func(:e1xb), Cvoid, (Ref{Float64}, Ref{Float64}), x, e1)
    e1[]
end

function _e1z(z::ComplexF64)
    ce1 = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE E1Z(Z,CE1)
    # double complex specfun_e1z(double complex z);
    ccall(f77func(:e1z), Cvoid, (Ref{ComplexF64}, Ref{ComplexF64}), z, ce1)
    ce1[]
end

function _eix(x::Float64)
    ei = Ref{Float64}(NaN)
    # SUBROUTINE EIX(X,EI)
    # double specfun_eix(double x);
    ccall(f77func(:eix), Cvoid, (Ref{Float64}, Ref{Float64}), x, ei)
    ei[]
end

function _eixz(z::ComplexF64)
    cei = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE EIXZ(Z,CEI)
    # double complex specfun_eixz(double complex z);
    ccall(f77func(:eixz), Cvoid, (Ref{ComplexF64}, Ref{ComplexF64}), z, cei)
    cei[]
end


"""Error Functions"""
#=Error function
- error
- cerror
- cerf
- cerzo
- cfc
- cfs
- fcs
- fcszo
- ffk
=#

function _erf(x::Float64)
    err = Ref{Float64}(NaN)
    # SUBROUTINE ERROR(X,ERR)
    ccall(f77func(:error), Cvoid, (Ref{Float64}, Ref{Float64},), x, err)
    err[]
end

function _erf(z::ComplexF64)
    cer = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE CERROR(Z,CER)
    ccall(f77func(:cerror), Cvoid, (Ref{ComplexF64}, Ref{ComplexF64}), z, cer)
    cer[]
end

function _cerf(z::ComplexF64)
    cer = Ref{ComplexF64}(NaN + NaN*im)
    cder = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE CERF(Z,CER,CDER)
    ccall(f77func(:cerf), Cvoid,
        (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
        z, cer, cder)
    cer[], cder[]
end

function _cerzo(zo::Vector{Complex{Float64}}, nt::Int)
    # SUBROUTINE CERZO(NT,ZO)
    ccall(f77func(:cerzo), Cvoid,
        (Ref{Int32}, Ref{ComplexF64}),
        nt, zo)
end

function _cfc(z::ComplexF64)
    zf = Ref{ComplexF64}(NaN + NaN*im)
    zd = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE CFC(Z,ZF,ZD)
    # void specfun_cfc(double complex z, double complex *zf, double complex *zd);
    ccall(f77func(:cfc), Cvoid,
        (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
        z, zf, zd)
    zf[], zd[]
end

function _cfs(z::ComplexF64)
    zf = Ref{ComplexF64}(NaN + NaN*im)
    zd = Ref{ComplexF64}(NaN + NaN*im)
    # SUBROUTINE CFS(Z,ZF,ZD)
    # void specfun_cfs(double complex z, double complex *zf, double complex *zd);
    ccall(f77func(:cfs), Cvoid,
        (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
        z, zf, zd)
    zf[], zd[]
end

function _fcs(z::Float64)
    c = Ref{Float64}(NaN)
    s = Ref{Float64}(NaN)
    # SUBROUTINE FCS(X,C,S)
    ccall(f77func(:fcs), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}),
        z, c, s)
    c[], s[]
end

function _fcszo(zo::Vector{Complex{Float64}}, kf::Int, nt::Int)
    # SUBROUTINE FCSZO(KF,NT,ZO)
    ccall(f77func(:fcszo), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{ComplexF64}),
        Int32(kf), Int32(nt), zo)
end

function _ffk(ks::Int, x::Float64)
    fr, fi, fm, fa, gr, gi, gm, ga = 
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN),
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN) 
    # SUBROUTINE FFK(KS,X,FR,FI,FM,FA,GR,GI,GM,GA)
    # void specfun_ffk(int ks, double x, 
    #   double *fr, double *fi, double *fm, double *fa,
    #   double *gr, double *gi, double *gm, double *ga);
    ccall(f77func(:ffk), Cvoid,
        (Ref{Int32}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        Int32(ks), x,
        fr, fi, fm, fa,
        gr, gi, gm, ga)
        
    fr[],fi[], fm[],fa[], gr[],gi[], gm[],ga[]
end


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
function _airyzo(
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
function _jynbh(
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
_jynbh(x, n::Int64, nmin::Int64, bj, by) =
    _jynbh(x, Int32(n), Int32(nmin), bj, by)

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
function _jyzo(n::Int64, nt::Int64,
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
    - _cgama

- [ ] _pbvv
    - [ ] _gamma2
    - [ ] _vvla
    - [ ] _vvsa

- [ ] _pbdv
    - [ ] _dvsa
    - [ ] _dvla

"""

function _pbwa(a::Float64, x::Float64)
    w1f, w2f, w1d, w2d =
        Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN), Ref{Float64}(NaN)
    # SUBROUTINE PBWA(A,X,W1F,W1D,W2F,W2D)
    ccall(f77func(:pbwa), Cvoid,
        (Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        a, x,
        w1f, w2f, w1d, w2d)
    w1f[], w2f[], w1d[], w2d[]
end

function _vvla(x::Float64, va::Float64)
    pv = Ref{Float64}(NaN)
    # SUBROUTINE VVLA(VA,X,PV)
    ccall(f77func(:vvla), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}),
        va, x, pv)
    pv[]
end

function _dvla(x::Float64, va::Float64)
    pd = Ref{Float64}(NaN)
    # SUBROUTINE DVLA(VA,X,PD)
    ccall(f77func(:dvla), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}),
        va, x, pd)
    pd[]
end

function _dvsa(x::Float64, va::Float64)
    pd = Ref{Float64}(NaN)
    # SUBROUTINE DVSA(VA,X,PD)
    ccall(f77func(:dvsa), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}),
        va, x, pd)
    pd[]
end

function _pbdv(dv::Vector{Float64}, dp::Vector{Float64}, x::Float64, v::Float64)
    pdf, pdd = Ref{Float64}(NaN), Ref{Float64}(NaN)
    # SUBROUTINE PBDV(V,X,DV,DP,PDF,PDD)
    ccall(f77func(:pbdv), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
         Ref{Float64}, Ref{Float64}),
        v, x, dv, dp,
        pdf, pdd)
    pdf[], pdd[]
end

function _vvsa(x::Float64, va::Float64)
    pv = Ref{Float64}(NaN)
    # SUBROUTINE VVSA(VA,X,PV)
    ccall(f77func(:vvsa), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}),
        va, x, pv)
    pv[]
end

function _pbvv(vv::Vector{Float64}, vp::Vector{Float64}, x::Float64, v::Float64)
    pvf, pvd = Ref{Float64}(NaN), Ref{Float64}(NaN)
    # SUBROUTINE PBVV(V,X,VV,VP,PVF,PVD)
    ccall(f77func(:pbvv), Cvoid,
        (Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
         Ref{Float64}, Ref{Float64}),
        v, x, vv, vp,
        pvf, pvd)
    pvf[], pvd[]
end


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


#=Legendre Functions
=#

"""
    SUBROUTINE LPMNS(M,N,X,  PM(0:N),PD(0:N))
    lpmns(int m, int n, T x,  T* pm, T* pd)

- Output: `pm, pd`
"""
function _lpmns(m::Int, n::Int, x::Float64, pm::Vector{Float64}, pd::Vector{Float64})
    ccall(f77func(:lpmns), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64},
         Ptr{Float64}, Ptr{Float64}),
        Int32(m), Int32(n), x,
        pm, pd)
end

"""
    SUBROUTINE LQMNS(M,N,X,  QM(0:N),QD(0:N))
    void lqmns(int m, int n, T x,  T *qm, T *qd)

- Output: `Qmn(x), Qmn'(x)`
"""
function _lqmns(m::Int, n::Int, x::Float64, qm::Vector{Float64}, qd::Vector{Float64})
    ccall(f77func(:lqmns), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64},
         Ptr{Float64}, Ptr{Float64}),
        Int32(m), Int32(n), x,
        qm, qd)
end

"""
    SUBROUTINE LPMV0(V,M,X,PMV)
    double lpmv0(double v, int m, double x)

- Output: `Pmv(x)`
"""
function _lpmv0(v::Float64, m::Int, x::Float64)
    pmv = Ref{Float64}(NaN)
    ccall(f77func(:lpmv0), Float64,
        (Ref{Float64}, Ref{Int32}, Ref{Float64},Ref{Float64}),
        v, Int32(m), x, pmv)
    pmv[]
end

"""
    SUBROUTINE LPMV(V,M,X, PMV)
    double lpmv(double x, int m, double v)

- Output: `Pmv(x)`
"""
function _lpmv(v::Float64, m::Int, x::Float64)
    pmv = Ref{Float64}(NaN)
    ccall(f77func(:lpmv), Float64,
        (Ref{Float64}, Ref{Int32}, Ref{Float64},Ref{Float64}),
        v, Int32(m), x, pmv)
    pmv[]
end


#= Spheroidal Wave Functions
- specfun_segv
- specfun_rswfp
- specfun_aswfa
- specfun_rswfo
=#

"""
Warp fortran `specfun.SEGV`.

- Input: `m,n,c`, `kd`
- Output: `cv`, `eg`
"""
function _segv(m::Int, n::Int, c::Float64, kd::Int, eg::Vector{Float64})
    cv = Ref{Float64}(NaN)
    @assert length(eg) >= 200
    # SUBROUTINE SEGV(M,N,C,KD,CV,EG)
    # void segv(int m, int n, T c, int kd, T *cv, T *eg)
    ccall(f77func(:segv), Cvoid,
        (Ref{Int32}, Ref{Int32},
         Ref{Float64}, Ref{Int32},
         Ref{Float64}, Ptr{Float64}),
         Int32(m), Int32(n),
         c, Int32(kd),
         cv, eg)
    cv[], eg
end

"""
    SUBROUTINE SDMN(M,N,C,CV,KD, DF(200))
    void sdmn(int m, int n, T c, T cv, int kd, T *df)

- Output: `df[]`
"""
function _sdmn(m::Int, n::Int, c::Float64, cv::Float64, kd::Int, df::Vector{Float64})
    ccall(f77func(:sdmn), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32},
         Ptr{Float64}),
        Int32(m), Int32(n), c, cv, Int32(kd),
        df)
    df
end

"""
    SUBROUTINE SCKB(M,N,C, DF(200),CK(200))
    void sckb(int m, int n, T c, T *df, T *ck)

- Input: `df[]`
- Output: `ck[]`
"""
function _sckb(m::Int, n::Int, c::Float64, df::Vector{Float64}, ck::Vector{Float64})
    ccall(f77func(:sckb), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64},
         Ptr{Float64}, Ptr{Float64}),
        Int32(m), Int32(n), c,
        df, ck)
    ck
end

"""
    SUBROUTINE ASWFA(M,N,C,X,KD,CV, S1F,S1D)
    void aswfa(int m, int n, T c, T x, int kd, T cv,  T *s1f, T *s1d)

- Output: `(s1f, s1d)`
"""
function _aswfa(m::Int, n::Int, c::Float64, x::Float64, kd::Int, cv::Float64) 
    s1f = Ref{Float64}(NaN)
    s1d = Ref{Float64}(NaN)
    ccall(f77func(:aswfa), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64},
         Ref{Float64}, Ref{Float64}),
        Int32(m), Int32(n), c, x, Int32(kd), cv,
        s1f, s1d)
    s1f[], s1d[]
end

"""
    SUBROUTINE SPHJ(N,X, NM,SJ(0:N),DJ(0:N))
    void sphj(T x, int n,  int *nm, T *sj, T *dj)

- Output: `(sj, dj, nm)`
"""
function _sphj(n::Int, x::Float64, sj::Vector{Float64}, dj::Vector{Float64})
    nm = Ref{Int32}(0)
    ccall(f77func(:sphj), Cvoid,
        (Ref{Int32}, Ref{Float64}, Ref{Int32},
         Ref{Float64}, Ref{Float64}),
        Int32(n), x, nm,
        sj, dj)
    sj, dj, Int(nm[])
end

"""
    SUBROUTINE RMN1(M,N,C,X,DF(200),KD, R1F,R1D)
    void rmn1(int m, int n, T c, T x, int kd, T *df, T *r1f, T *r1d)

- Output: `(r1f, r1d)`
"""
function _rmn1(m::Int, n::Int, c::Float64, x::Float64, kd::Int, df::Vector{Float64})
    r1f, r1d  = Ref{Float64}(NaN), Ref{Float64}(NaN)
    ccall(f77func(:rmn1), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Int32},
         Ref{Float64}, Ref{Float64}),
        Int(m), Int(n), c, x, df, Int(kd),
        r1f, r1d)
    r1f[], r1d[]
end

"""
    SUBROUTINE SPHY(N,X,NM, SY(0:N),DY(0:N))
    void sphy(T x, int n, int *nm, T *sy, T *dy)

- Output: `(sy, dy, nm)`
"""
function _sphy(n::Int, x::Float64, sy::Vector{Float64}, dy::Vector{Float64})
    nm = Ref{Int32}(0)
    ccall(f77func(:sphy), Cvoid,
        (Ref{Int32}, Ref{Float64}, Ref{Int32},
         Ptr{Float64}, Ptr{Float64}),
        Int32(n), x, nm,
        sy, dy)
    sy, dy, Int(nm[])
end

"""
    SUBROUTINE RMN2L(M,N,C,X, DF(200), KD,R2F,R2D,ID)
    void rmn2l(int m, int n, T c, T x, int Kd, T *Df,  T *R2f, T *R2d, int *Id)

- Output: `(r2f, r2d, id)`
"""
function _rmn2l(m::Int, n::Int, c::Float64, x::Float64, kd::Int, df::Vector{Float64})
    r2f = Ref{Float64}(NaN)
    r2d = Ref{Float64}(NaN)
    id = Ref{Int32}(0)
    ccall(f77func(:rmn2l), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
         Ptr{Float64}, Ref{Int32},
         Ref{Float64}, Ref{Float64}, Ref{Int32}),
        Int32(m), Int32(n), c, x,
        df, Int32(kd),
        r2f, r2d, id)
    r2f[], r2d[], Int(id[])
end
