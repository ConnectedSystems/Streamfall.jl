push!(LOAD_PATH,"../src/")

# using Pkg

# pkg"activate .."

using Documenter, Streamfall


makedocs(sitename="Streamfall Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "index.md",
        "Streamfall.md",
        "metrics.md",
        "Examples" => [
            "Calibration" => [
                "calibration_setup.md",
                "calibration.md",
            ],
            "simple_showcase.md",
            "multisystem_showcase.md",
        ]
        
    ]
)

deploydocs(
    repo = "github.com/ConnectedSystems/Streamfall.jl.git",
    devbranch = "main"
)