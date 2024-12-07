# PureSpecial

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://inkydragon.github.io/PureSpecial.jl/dev/)
[![Build Status](https://github.com/inkydragon/PureSpecial.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/inkydragon/PureSpecial.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/inkydragon/PureSpecial.jl/graph/badge.svg?token=krUDB5Fepa)](https://codecov.io/gh/inkydragon/PureSpecial.jl)

> Rewriting scipy's fortran dependencies in julia.  
> Currently focusing on `scipy.special`


## Test

Activate test project:

```sh
julia --project=test
```

Run test and gen coverage:

```jl
using Pkg; using LocalCoverage;
ENV["GITHUB_WORKSPACE"] = pwd()
Pkg.add(url=".");  html_coverage(generate_coverage("PureSpecial"; run_test=true); dir = "../cov")
# Pkg.add(PackageSpec(path=".", rev="dev")); html_coverage(generate_coverage("PureSpecial"; run_test=true); dir = "../cov")
```


## License

```jl
# SPDX-License-Identifier: MIT
```
