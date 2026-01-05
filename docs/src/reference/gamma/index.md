# Gamma Functions

```@docs
PureSpecial.ASF.GammaFunctionsDoc
```

## Factorial Function

| Function              | Description               |
|:----------------------|:--------------------------|
| `Base.factorial`      | factorial function    |
| [`factorial2`](@ref)  | double factorial function |
| [`factorialk`](@ref)  | multifactorial function   |
| `Base.binomial`       | binomial coefficient  |

```@docs
PureSpecial.ASF.factorial
PureSpecial.ASF.doublefactorial
PureSpecial.ASF.multifactorial
PureSpecial.ASF.superfactorial
PureSpecial.ASF.hyperfactorial
```

### Binomial Coefficient

```@docs
PureSpecial.ASF.binomial
PureSpecial.ASF.multinomial
```


## Gamma Function

| Function              | Description               |
|:----------------------|:--------------------------|
| [`gamma`](@ref)       | gamma function            |
| [`loggamma`](@ref)    | log gamma function        |
| [`logabsgamma`](@ref) | log abs gamma function    |
| [`digamma`](@ref)     | digamma function          |
| [`trigamma`](@ref)    | trigamma function         |
| [`polygamma`](@ref)   | polygamma function        |

```@docs
PureSpecial.ASF.gamma
PureSpecial.ASF.digamma
PureSpecial.ASF.polygamma
```

## Incomplete Gamma Function

## Pochhammer Function

## Beta Function

```@docs
PureSpecial.ASF.beta
```

## Autodoc

```@autodocs
Modules = [PureSpecial]
Pages   = ["gamma.jl"]
```
