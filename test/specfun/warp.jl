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
