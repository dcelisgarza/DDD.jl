# IO
## Input
```@docs
load(filename::AbstractString)

loadDislocationLoop(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)

loadMaterialP(dict::Dict{T1, T2}) where {T1, T2}

loadIntegrationP(dict::Dict{T1, T2}) where {T1, T2}

loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}

loadDislocationP(dict::Dict{T1, T2}) where {T1, T2}

loadParams(
    fileDislocationP::AbstractString,
    fileMaterialP::AbstractString,
    fileIntegrationP::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)

loadNetwork(fileDislocationNetwork::AbstractString)
```

## Output

```@docs
save(filename::AbstractString, args...; mode::AbstractString = "w")
```
