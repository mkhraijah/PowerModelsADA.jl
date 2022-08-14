push!(LOAD_PATH,"../src/")
using Documenter, PMADA

makedocs(
    modules = [PMADA],
    sitename = "PMADA.jl",
    authors = "Mohannad Alkhraijah",
    format = Documenter.HTML(
        analytics = "",
        mathengine = Documenter.MathJax(),
        collapselevel=1,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Quick Start Guide" => "quickguide.md",
            "Distributed Algorithms" => [
                "ADMM" => "admm.md"
            ]
        ]
    ]

)

deploydocs(
    repo = "github.com/mkhraijah/PMADA.jl"
)
