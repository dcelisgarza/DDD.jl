# IO
## Input
```@docs
loadJSON(filename::AbstractString)

loadDislocationLoopJSON(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)

loadMaterialParametersJSON(dict::Dict{T1, T2}) where {T1, T2}

loadIntegrationParametersJSON(dict::Dict{T1, T2}) where {T1, T2}

loadSlipSystemJSON(dict::Dict{T1, T2}) where {T1, T2}

loadDislocationParametersJSON(dict::Dict{T1, T2}) where {T1, T2}

loadParametersJSON(
    fileDislocationParameters::AbstractString,
    fileMaterialParameters::AbstractString,
    fileIntegrationParameters::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)

loadNetworkJSON(fileDislocationNetwork::AbstractString)
```

## Output

```@docs
saveJSON(filename::AbstractString, args...; mode::AbstractString = "w")
```
