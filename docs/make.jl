using PureSpecial
using Documenter

DocMeta.setdocmeta!(PureSpecial, :DocTestSetup, :(using PureSpecial); recursive=true)


pages = Any[
    "Home" => "index.md",
    "Docs" => Any[
        "ref/struve.md"
    ],
    "specfun.f" => Any[
        "specfun/index.md",
        "specfun/book-chapters.md",
        "specfun/book-index.md",
        "specfun/autodocs.md",
    ],
    "cephes/index.md",
    "faddeeva/index.md",
    "impls.md",
    "scipy.special.md",
    "special-functions.md",
]

makedocs(;
    modules=[PureSpecial, PureSpecial.Specfun],
    repo=Remotes.GitHub("inkydragon", "PureSpecial.jl"),
    authors="Chengyu HAN <cyhan.dev@outlook.com> and contributors",
    sitename="PureSpecial.jl",
    format=Documenter.HTML(;
        repolink="https://inkydragon.github.io/PureSpecial.jl",
        canonical="https://inkydragon.github.io/PureSpecial.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/inkydragon/PureSpecial.jl",
    devbranch="main",
)
