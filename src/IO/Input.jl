"""
```
load(filename::AbstractString)
```
Wrapper for `JSON.parsefile(filename)`.
"""
function load(filename::AbstractString)
    dict = JSON.parsefile(filename)
    return dict
end

"""
```
function loadDislocationLoop(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)
```
Loads initial dislocation structure out of a dictionary loaded from a JSON file. Returns a variable of type [`DislocationLoop`](@ref).
"""
function loadDislocationLoop(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)

    dlnTypes = makeTypeDict(AbstractDlnStr)
    distributions = makeTypeDict(AbstractDistribution)

    slipPlane = slipSystem.slipPlane
    bVec = slipSystem.bVec

    range = zeros(Float64, 2, 3)
    for i = 1:3
        range[:, i] = convert.(Int64, dict["range"][i])
    end

    dislocationLoop = DislocationLoop(
        loopType = dlnTypes[dict["loopType"]],
        numSides = convert(Int64, dict["numSides"]),
        nodeSide = convert(Int64, dict["nodeSide"]),
        numLoops = convert(Int64, dict["numLoops"]),
        segLen = convert.(Float64, dict["segLen"]),
        slipSystem = convert.(Int64, dict["slipSystem"]),
        _slipPlane = convert.(Float64, slipPlane[dict["slipSystem"], 1:3]),
        _bVec = convert.(Float64, bVec[dict["slipSystem"], 1:3]),
        label = nodeType.(dict["label"]),
        buffer = convert(Float64, dict["buffer"]),
        range = range,
        dist = distributions[dict["dist"]],
    )

    return dislocationLoop
end

"""
```
loadMaterialP(dict::Dict{T1, T2}) where {T1, T2}
```
Loads material parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`MaterialP`](@ref).
"""
function loadMaterialP(dict::Dict{T1, T2}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    materialP = MaterialP(
        μ = convert(Float64, dict["μ"]),
        μMag = convert(Float64, dict["μMag"]),
        ν = convert(Float64, dict["ν"]),
        E = convert(Float64, dict["E"]),
        crystalStruct = crystalStruct[dict["crystalStruct"]],
    )

    return materialP
end

"""
```
loadIntegrationP(dict::Dict{T1, T2}) where {T1, T2}
```
Loads integration parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`IntegrationP`](@ref).
"""
function loadIntegrationP(dict::Dict{T1, T2}) where {T1, T2}

    integDict = makeTypeDict(AbstractIntegrator)

    integrationP = IntegrationP(
        dt = convert(Float64, dict["dt"]),
        tmin = convert(Float64, dict["tmin"]),
        tmax = convert(Float64, dict["tmax"]),
        method = integDict[dict["method"]],
        abstol = convert(Float64, dict["abstol"]),
        reltol = convert(Float64, dict["reltol"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int64, dict["step"]),
    )

    return integrationP
end

"""
```
loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}
```
Loads slip systems out of a dictionary loaded from a JSON file. Returns a variable of type [`SlipSystem`](@ref).
"""
function loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    slipSystem = SlipSystem(
        crystalStruct = crystalStruct[dict["crystalStruct"]],
        slipPlane = convert.(
            Float64,
            [dict["slipPlane"][1] dict["slipPlane"][2] dict["slipPlane"][3]],
        ),
        bVec = convert.(
            Float64,
            [dict["bVec"][1] dict["bVec"][2] dict["bVec"][3]],
        ),
    )

    return slipSystem
end

"""
```
loadDislocationP(dict::Dict{T1, T2}) where {T1, T2}
```
Loads dislocation parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`DislocationP`](@ref).
"""
function loadDislocationP(dict::Dict{T1, T2}) where {T1, T2}

    mobDict = makeTypeDict(AbstractMobility)

    dislocationP = DislocationP(
        coreRad = convert(Float64, dict["coreRad"]),
        coreRadMag = convert(Float64, dict["coreRadMag"]),
        minSegLen = convert(Float64, dict["minSegLen"]),
        maxSegLen = convert(Float64, dict["maxSegLen"]),
        minArea = convert(Float64, dict["minArea"]),
        maxArea = convert(Float64, dict["maxArea"]),
        maxConnect = convert(Int64, dict["maxConnect"]),
        remesh = dict["remesh"],
        collision = dict["collision"],
        separation = dict["separation"],
        virtualRemesh = dict["virtualRemesh"],
        edgeDrag = convert(Float64, dict["edgeDrag"]),
        screwDrag = convert(Float64, dict["screwDrag"]),
        climbDrag = convert(Float64, dict["climbDrag"]),
        lineDrag = convert(Float64, dict["lineDrag"]),
        mobility = mobDict[dict["mobility"]],
    )

    return dislocationP
end

"""
```
loadParams(
    fileDislocationP::AbstractString,
    fileMaterialP::AbstractString,
    fileIntegrationP::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)
```
Loads simulation parameters out of a dictionary loaded from a JSON file. Returns a tuple of variable types ([`DislocationP`](@ref), [`MaterialP`](@ref), [`IntegrationP`](@ref), [`SlipSystem`](@ref), [`DislocationLoop`](@ref)) or vectors of those types.
"""
function loadParams(
    fileDislocationP::AbstractString,
    fileMaterialP::AbstractString,
    fileIntegrationP::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)
    # We use JSON arrays because it lets us dump a variable number of args into a single JSON file. To keep things gonsistent we use them always. Hence the indices here.
    dictDislocationP = load(fileDislocationP)
    dislocationP = loadDislocationP(dictDislocationP[1])

    dictMaterialP = load(fileMaterialP)
    materialP = loadMaterialP(dictMaterialP[1])

    dictIntegrationP = load(fileIntegrationP)
    integrationP = loadIntegrationP(dictIntegrationP[1])

    dictSlipSystem = load(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem[1])
    # There can be multiple dislocations per simulation parameters.
    dictDislocationLoop = load(fileDislocationLoop)
    dislocationLoop = zeros(DislocationLoop, length(dictDislocationLoop))
    for i in eachindex(dislocationLoop)
        dislocationLoop[i] =
            loadDislocationLoop(dictDislocationLoop[i], slipSystems)
    end

    return dislocationP, materialP, integrationP, slipSystems, dislocationLoop
end # function

"""
```
loadNetwork(fileDislocationNetwork::AbstractString)
```
Loads a dislocation network from a JSON file.
"""
function loadNetwork(fileDislocationNetwork::AbstractString)
    dict = load(fileDislocationNetwork)[1]

    lenLinks = length(dict["links"][1])
    lenCoord = length(dict["coord"][1])
    maxConnect = convert.(Int64, dict["maxConnect"])

    links = zeros(Int64, lenLinks, 2)
    slipPlane = zeros(Float64, lenLinks, 3)
    bVec = zeros(Float64, lenLinks, 3)
    coord = zeros(Float64, lenCoord, 3)
    connectivity = zeros(Int64, lenLinks, 2 * maxConnect + 1)
    linksConnect = zeros(Int64, lenLinks, 2)
    segIdx = zeros(Int64, lenLinks, 3)

    for i = 1:2
        links[:, i] = convert.(Int64, dict["links"][i])
        linksConnect[:, i] = convert.(Int64, dict["linksConnect"][i])
    end

    @inbounds @simd for i = 1:3
        slipPlane[:, i] = convert.(Float64, dict["slipPlane"][i])
        bVec[:, i] = convert.(Float64, dict["bVec"][i])
        coord[:, i] = convert.(Float64, dict["coord"][i])
        segIdx[:, i] = convert.(Float64, dict["segIdx"][i])
    end

    @inbounds @simd for i = 1:9
        connectivity[:, i] = convert.(Int64, dict["connectivity"][i])
    end

    dislocationNetwork = DislocationNetwork(
        links = links,
        slipPlane = slipPlane,
        bVec = bVec,
        coord = coord,
        label = nodeType.(dict["label"]),
        numNode = convert.(Int64, dict["numNode"]),
        numSeg = convert.(Int64, dict["numSeg"]),
        maxConnect = maxConnect,
    )
    dislocationNetwork.linksConnect = linksConnect
    dislocationNetwork.connectivity = connectivity
    dislocationNetwork.segIdx = segIdx

    return dislocationNetwork
end
