push!(LOAD_PATH,"../src/")
using Documenter
using PMADA

makedocs(
    sitename = "PMADA.jl",
    format = Documenter.HTML(
        analytics = "",
        mathengine = Documenter.MathJax(),
        prettyurls=false,
        collapselevel=1,
    ),
    modules = [PMADA],
    authors = "Mohannad Alkhraijah",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mkhraijah/PMADA.jl"
)
