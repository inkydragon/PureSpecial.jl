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
- fresnelf, fresnelg

## Voigt Function
- voigtu, voigtv

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

