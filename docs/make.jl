using SymbolicCRN
using Documenter

DocMeta.setdocmeta!(SymbolicCRN, :DocTestSetup, :(using SymbolicCRN); recursive=true)

makedocs(;
    modules=[SymbolicCRN],
    authors="Laura Brustenga i Moncusí <brust@math.ku.dk> and contributors",
    repo="https://github.com/LauraBMo/SymbolicCRN.jl/blob/{commit}{path}#{line}",
    sitename="SymbolicCRN.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LauraBMo.github.io/SymbolicCRN.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/SymbolicCRN.jl",
)
