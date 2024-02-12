using Scipy4j
using Documenter

DocMeta.setdocmeta!(Scipy4j, :DocTestSetup, :(using Scipy4j); recursive=true)

makedocs(;
    modules=[Scipy4j],
    authors="Chengyu HAN <git@wo-class.cn> and contributors",
    sitename="Scipy4j.jl",
    format=Documenter.HTML(;
        canonical="https://inkydragon.github.io/Scipy4j.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/inkydragon/Scipy4j.jl",
    devbranch="main",
)
