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
        "network.md",
        "nodes.md",
        "use_methods.md",
        "Examples" => [
            "Calibration" => [
                "metrics.md",
                "calibration_setup.md",
                "calibration.md",
            ],
            "simple_showcase.md",
            # "multisystem_showcase.md",
        ]
        
    ]
)

deploydocs(
    repo = "github.com/ConnectedSystems/Streamfall.jl.git",
    devbranch = "main",
    target="build",
    deps=nothing,
    make=nothing
)
