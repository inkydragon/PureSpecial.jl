# SPDX-License-Identifier: MIT
"""
# Exponential and Trigonometric Integrals

## Exponential Integral

- [`expint1`](@ref)
- [`expint`](@ref), [`expintx`](@ref)
- [`expinti`](@ref)

### Logarithmic integral

- [`logint`](@ref)

## Trigonometric Integral

- [`sinint`](@ref), [`cosint`](@ref)
- [`sinhint`](@ref), [`coshint`](@ref)

# Reference

- [DLMF: Chapter 6 Exponential, Logarithmic, Sine, and Cosine Integrals](https://dlmf.nist.gov/6)
"""
const ExponentialIntegralsDoc = nothing


"""
    expint1(z)

(principal value of) exponential integral, `E₁(z)`.

Also known as:
`exp1(z)`, `expint(z)`, `expint_E1`, `ExpIntegralE[1,z]`
"""
function expint1 end 

"""
    expint(ν, z)

Generalized exponential integral, `En(z)`, `Eν(z)`.

Also known as:
`expint(n,z)`, `expint_En`, `ExpIntegralE[n,z]`
"""
function expint end 

"""
    expintx(ν, z)

Scaled (generalized) exponential integral, `eᶻEν(z)`.
"""
function expintx end 

"""
    expinti(z)

Exponential integral, `Ei(z)`.

Also known as:
`expi(z)`, `expint_Ei`, `ExpIntegralEi[z]`
"""
function expinti end 

"""
    logint(z)

Logarithmic integral, `Li(z)`.

Also known as:
`logint(z)`, `LogIntegral[z]`
"""
function logint end 

"""
    sinint(z)

Sine integral, `Si(z)`.

Also known as:
`Si`, `sici(z)`, `SinIntegral[z]`
"""
function sinint end 

"""
    cosint(z)

Cosine integral, `Ci(z)`.

Also known as:
`Ci`, `sici(z)`, `CosIntegral[z]`
"""
function cosint end 

"""
    sinhint(z)

Hyperbolic sine integral, `Shi(z)`.

Also known as:
`Shi`, `shichi(z)`, `SinhIntegral[z]`
"""
function sinhint end 

"""
    coshint(z)

Hyperbolic cosine integral, `Chi(z)`.

Also known as:
`Chi`, `shichi(z)`, `CoshIntegral[z]`
"""
function coshint end 
