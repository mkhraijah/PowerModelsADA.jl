push!(LOAD_PATH,"../src/")

using PMADA
using Documenter

makedocs(sitename="Documentation")

deploydocs(
    repo = "github.gatech.edu/malkhraijah3/PMADA.jl.git",
)
