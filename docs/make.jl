using Documenter, DDD
using DDD: removeNode!, removeConnection!, removeLink!

makedocs(;
    modules = [DDD],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    repo = "https://github.com/dcelisgarza/DDD.jl/blob/{commit}{path}#L{line}",
    sitename = "DDD.jl",
    authors = "Daniel Celis Garza",
    # assets=String[],
)

deploydocs(; repo = "github.com/dcelisgarza/DDD.jl.git")
