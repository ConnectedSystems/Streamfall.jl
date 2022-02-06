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
        "primer.md",
        "network.md",
        "Nodes" => [
            "Node.md",
            "IHACRES.md",
            "HyMod.md",
            "GR4J.md",
            "SYMHYD.md",
            "Dam.md"
        ],
        "use_methods.md",
        "metrics.md",
        "Examples" => [
            "Calibration" => [
                "calibration_setup.md",
                "calibration.md",
            ],
            "simple_showcase.md",
            "model_comparison.md",
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
