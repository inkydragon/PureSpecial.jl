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
