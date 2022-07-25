push!(LOAD_PATH,"../src/")

using PMADA
using Documenter

makedocs(
    modules = [PMADA],
    format = Documenter.HTML(
        analytics = "",
        mathengine = Documenter.MathJax(),
        prettyurls=false,
        collapselevel=1,
    ),
    strict=false,
    sitename = "PMADA.jl",
    authors = "Mohannad Alkhraijah",
    pages = ["Introduction" => "index.md"]
)

deploydocs(
    repo = "github.com/mkhraijah/PMADA.jl.git",
)
