using PureSpecial
using Documenter

DocMeta.setdocmeta!(PureSpecial, :DocTestSetup, :(using PureSpecial); recursive=true)


pages = Any[
    "Home" => "index.md",
    "Special functions" => "special.md",
]

makedocs(;
    modules=[PureSpecial],
    repo=Remotes.GitHub("inkydragon", "PureSpecial.jl"),
    authors="Chengyu HAN <cyhan.dev@outlook.com> and contributors",
    sitename="PureSpecial.jl",
    format=Documenter.HTML(;
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
