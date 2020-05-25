using Documenter, DDD

import DDD:
    AbstractMesh,
    RegularCuboidMesh,
    DislocationFEMCorrective,
    calcPKForce,
    makeInstanceDict,
    translateEnum,
    subTypeTree,
    inclusiveComparison,
    dimDot,
    dimNorm,
    AbstractCrystalStruct,
    MaterialP,
    AbstractIntegrator,
    IntegrationP,
    AbstractShapeFunction
makedocs(;
    modules = [DDD],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Types" => "types.md",
        "IO" => "io.md",
        "Post Processing" => "postProcessing.md",
        "Functions" => "functions.md",
        "Discrete Dislocation Dynamics" => "theory.md",
        "Motivation" => "motivation.md",
        "Index" => "idx.md",
    ],
    repo = "https://github.com/dcelisgarza/DDD.jl/blob/{commit}{path}#L{line}",
    sitename = "DDD.jl",
    authors = "Daniel Celis Garza",
    # assets=String[],
)

deploydocs(; repo = "github.com/dcelisgarza/DDD.jl.git")
