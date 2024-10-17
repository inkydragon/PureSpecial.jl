using Scipy4j
using Documenter

DocMeta.setdocmeta!(Scipy4j, :DocTestSetup, :(using Scipy4j); recursive=true)


pages = Any[
    "Home" => "index.md",
    "Special functions" => "special.md",
]

makedocs(;
    modules=[Scipy4j],
    repo=Remotes.GitHub("inkydragon", "Scipy4j.jl"),
    authors="Chengyu HAN <cyhan.dev@outlook.com> and contributors",
    sitename="Scipy4j.jl",
    format=Documenter.HTML(;
        canonical="https://inkydragon.github.io/Scipy4j.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/inkydragon/Scipy4j.jl",
    devbranch="main",
)
