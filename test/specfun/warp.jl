# SPDX-License-Identifier: MIT
using Libdl

"""
## Compile fortran lib
```sh
gfortran -I. -shared -fPIC -g -O3 specfun.f -o libspecfun.so
```
"""
# where to find: `specfun.so`
push!(DL_LOAD_PATH, "/mnt/d/jl/JuliaMath/Scipy4j.jl")
push!(DL_LOAD_PATH, pwd())
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
        (Ref{Int64},
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
