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
            "Problem Specification and Algorithm Selection" => "specification.md" , 
            "Distributed Algorithms" => [
                "ADMM" => "admm.md",
                "ATC" => "atc.md",
                "APP" => "app.md",
                "ALADEN" => "aladin.md"
            ]
        ],
        "Tutorials" => [
            "Using PowerModelsADA" => "tutorial.md",
        "User-defined Algorithm" => "newalgorithm.md"
        ],
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/mkhraijah/PowerModelsADA.jl"
)
