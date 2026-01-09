# PureSpecial

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://inkydragon.github.io/PureSpecial.jl/dev/)
[![Build Status](https://github.com/inkydragon/PureSpecial.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/inkydragon/PureSpecial.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/inkydragon/PureSpecial.jl/graph/badge.svg?token=krUDB5Fepa)](https://codecov.io/gh/inkydragon/PureSpecial.jl)

> Rewriting [scipy/xsf](https://github.com/scipy/xsf) in julia.  


## Dev

### Build docs

```sh
# The following command will init docs project in the `docs/` directory.
#   You only need to run this line once.
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate();"
julia --project=docs docs/make.jl
# html files located in `docs/build/`
```

### Test

> You need a fortran compiler to run the test.

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
