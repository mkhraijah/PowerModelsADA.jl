push!(LOAD_PATH,"../src/")
using Documenter, PowerModelsADA

makedocs(
    modules = [PowerModelsADA],
    sitename = "PowerModelsADA.jl",
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
            "Data Structure" => "data_structure.md",
            "Technical Specifications" => "specification.md" , 
            "Distributed Algorithms" => [
                "ADMM" => "admm.md",
                "ATC" => "atc.md",
                "APP" => "app.md",
                "ALADIN" => "aladin.md",
                "Adaptive ADMM" => "adaptive_admm.md"
            ]
        ],
        "Tutorials" => [
            "Using PowerModelsADA" => "tutorial.md",
        "User-defined Algorithm" => "newalgorithm.md"
        ],
        "Library" => "library.md",
        "Comparison Results" => "comparison.md"

    ]
)

deploydocs(
    repo = "github.com/mkhraijah/PowerModelsADA.jl"
)
