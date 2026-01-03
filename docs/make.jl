using PureSpecial
using Documenter

DocMeta.setdocmeta!(PureSpecial, :DocTestSetup, :(using PureSpecial); recursive=true)


reference_pages = Any[
    "Gamma Functions" => Any[
        "reference/gamma/index.md"
        "reference/gamma/gamma.md"
    ],
    "reference/exp-integral/index.md",

    # TODO: plot is slow
    # "Struve Functions" => Any[
    #     "reference/struve.md",
    # ],
    "reference/autodocs.md",
] # reference_pages

pages = Any[
    "Home" => "index.md",
    "special-functions.md",
    "Reference" => reference_pages,

    "Dev doc" => Any[
        "impls.md",
        "scipy.special.md",
        "specfun.f" => Any[
            "specfun/index.md",
            "specfun/book-chapters.md",
            "specfun/book-index.md",
            "specfun/autodocs.md",
        ],
        "cephes/index.md",
        "faddeeva/index.md",
    ],
]

makedocs(;
    modules=[PureSpecial, PureSpecial.Specfun],
    repo=Remotes.GitHub("inkydragon", "PureSpecial.jl"),
    authors="Chengyu HAN <cyhan.dev@outlook.com> and contributors",
    sitename="PureSpecial.jl",
    format=Documenter.HTML(;
        # canonical="https://inkydragon.github.io/PureSpecial.jl",
        canonical="https://cyhan.dev/PureSpecial.jl",
        edit_link="main",
        assets=String[],
    ),
    warnonly=true,
    checkdocs=:exports,
    pages=pages,
)

deploydocs(;
    repo="github.com/inkydragon/PureSpecial.jl",
    devbranch="main",
    push_preview=true,
)
