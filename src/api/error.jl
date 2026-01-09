# SPDX-License-Identifier: MIT
"""
# Error Functions

## Error Function
- [`erf`](@ref), [`erfc`](@ref), [`erfi`](@ref)
- [`faddeeva(z)`](@ref)

## Dawson Integral
- [`dawson(z)`](@ref)

## Fresnel Integral
- [`fresnelc(z)`](@ref), [`fresnels(z)`](@ref)
- [`fresnelf(z)`](@ref), [`fresnelg(z)`](@ref)
- `ℱ(z)`, `G(z)`

## Voigt Function
- [`voigtu(x,t)`](@ref), [`voigtv(x,t)`](@ref)

# Reference
- [DLMF: Chapter 7 Error Functions, Dawson’s and Fresnel Integrals](https://dlmf.nist.gov/7)
"""
const ErrorFunctionsDoc = nothing


"""
    erf(x::Real)
    erf(z::Complex)

Error function.
"""
function erf end

"""
    erfc(x::Real)
    erfc(z::Complex)

Complementary error function.
"""
function erfc end

"""
    erfi(x::Real)
    erfi(z::Complex)

Imaginary error function.
"""
function erfi end

"""
    faddeeva(x::Real)
    faddeeva(z::Complex)

Faddeeva function, `w(z)`.
"""
function faddeeva end

"""
    dawson(x::Real)
    dawson(z::Complex)

Dawson integral.
"""
function dawson end

"""
    fresnelc(x::Real)
    fresnelc(z::Complex)

Fresnel integral, `C(z)`.
"""
function fresnelc end

"""
    fresnels(x::Real)
    fresnels(z::Complex)

Fresnel integral, `S(z)`.
"""
function fresnels end

"""
    fresnelf(x::Real)
    fresnelf(z::Complex)

Auxiliary function for Fresnel integrals, `f(z)`.
"""
function fresnelf end

"""
    fresnelg(x::Real)
    fresnelg(z::Complex)

Auxiliary function for Fresnel integrals, `g(z)`.
"""
function fresnelg end

"""
    voigtu(x::Real, t::Real)

Voigt function, `U(x,t)`.
"""
function voigtu end

"""
    voigtv(x::Real, t::Real)

Voigt function, `V(x,t)`.
"""
function voigtv end
