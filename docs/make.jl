using Documenter, DDD

makedocs(;
    modules = [DDD],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Types" => "Types.md",
        "Constructors" => "Constructors.md"
        # "Dislocations" => "Dislocations.md",
        # "IO" => "io.md",
        # "Post Processing" => "postProcessing.md",
        # "Functions" => "functions.md",
        # "Discrete Dislocation Dynamics" => "theory.md",
        # "Motivation" => "motivation.md",
        # "Theory" => "theory.md",
        # "Index" => "idx.md",
    ],
    repo = "https://github.com/dcelisgarza/DDD.jl/blob/{commit}{path}#L{line}",
    sitename = "DDD.jl",
    authors = "Daniel Celis Garza",
    # assets=String[],
)

deploydocs(; repo = "github.com/dcelisgarza/DDD.jl.git")
