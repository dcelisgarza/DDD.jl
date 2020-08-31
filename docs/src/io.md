# IO
## Input
```@docs
load(filename::AbstractString)

loadDislocationLoop(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)

loadMaterialParameters(dict::Dict{T1, T2}) where {T1, T2}

loadIntegrationParameters(dict::Dict{T1, T2}) where {T1, T2}

loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}

loadDislocationParameters(dict::Dict{T1, T2}) where {T1, T2}

loadParams(
    fileDislocationParameters::AbstractString,
    fileMaterialParameters::AbstractString,
    fileIntegrationParameters::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)

loadNetwork(fileDislocationNetwork::AbstractString)
```

## Output

```@docs
save(filename::AbstractString, args...; mode::AbstractString = "w")
```
