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
function loadDislocationLoop(dict::Dict{T1, T2} where {T1, T2}, slipSystem::SlipSystem)

    dlnTypes = makeTypeDict(AbstractDlnStr)
    distributions = makeTypeDict(AbstractDistribution)

    slipPlane = slipSystem.slipPlane
    bVec = slipSystem.bVec

    range = zeros(3, 2)
    for i in 1:2
        range[:, i] = convert.(Int, dict["range"][i])
    end

    dislocationLoop = DislocationLoop(;
        loopType = dlnTypes[dict["loopType"]],
        numSides = convert(Int, dict["numSides"]),
        nodeSide = convert(Int, dict["nodeSide"]),
        numLoops = convert(Int, dict["numLoops"]),
        segLen = convert.(Float64, dict["segLen"]),
        slipSystem = convert.(Int, dict["slipSystem"]),
        _slipPlane = convert.(Float64, slipPlane[:, dict["slipSystem"]]),
        _bVec = convert.(Float64, bVec[:, dict["slipSystem"]]),
        label = nodeType.(dict["label"]),
        buffer = convert(Float64, dict["buffer"]),
        range = range,
        dist = distributions[dict["dist"]],
    )

    return dislocationLoop
end

"""
```
loadMaterialParameters(dict::Dict{T1, T2}) where {T1, T2}
```
Loads material parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`MaterialParameters`](@ref).
"""
function loadMaterialParameters(dict::Dict{T1, T2}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    MaterialParams = MaterialParameters(
        μ = convert(Float64, dict["μ"]),
        μMag = convert(Float64, dict["μMag"]),
        ν = convert(Float64, dict["ν"]),
        E = convert(Float64, dict["E"]),
        crystalStruct = crystalStruct[dict["crystalStruct"]],
        σPN = convert(Float64, dict["σPN"]),
    )

    return MaterialParams
end

"""
```
loadIntegrationParameters(dict::Dict{T1, T2}) where {T1, T2}
```
Loads integration parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`IntegrationParameters`](@ref).
"""
function loadIntegrationParameters(dict::Dict{T1, T2}) where {T1, T2}

    integDict = makeTypeDict(AbstractIntegrator)

    IntegrationParams = IntegrationParameters(;
        method = integDict[dict["method"]],
        tmin = convert(Float64, dict["tmin"]),
        tmax = convert(Float64, dict["tmax"]),
        dtmin = convert(Float64, dict["dtmin"]),
        dtmax = convert(Float64, dict["dtmax"]),
        abstol = convert(Float64, dict["abstol"]),
        reltol = convert(Float64, dict["reltol"]),
        maxchange = convert(Float64, dict["maxchange"]),
        exponent = convert(Float64, dict["exponent"]),
        maxiter = convert(Int, dict["exponent"]),
    )

    return IntegrationParams
end

"""
```
loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}
```
Loads slip systems out of a dictionary loaded from a JSON file. Returns a variable of type [`SlipSystem`](@ref).
"""
function loadSlipSystem(dict::Dict{T1, T2}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    lenSlipSys = length(dict["slipPlane"])
    slipPlane = zeros(3, lenSlipSys)
    bVec = zeros(3, lenSlipSys)

    for i in 1:lenSlipSys
        slipPlane[:, i] = convert.(Float64, dict["slipPlane"][i])
        bVec[:, i] = convert.(Float64, dict["bVec"][i])
    end

    slipSystem = SlipSystem(;
        crystalStruct = crystalStruct[dict["crystalStruct"]],
        slipPlane = slipPlane,
        bVec = bVec,
    )

    return slipSystem
end

"""
```
loadDislocationParameters(dict::Dict{T1, T2}) where {T1, T2}
```
Loads dislocation parameters out of a dictionary loaded from a JSON file. Returns a variable of type [`DislocationParameters`](@ref).
"""
function loadDislocationParameters(dict::Dict{T1, T2}) where {T1, T2}

    mobDict = makeTypeDict(AbstractMobility)

    DislocationParams = DislocationParameters(;
        coreRad = convert(Float64, dict["coreRad"]),
        coreRadMag = convert(Float64, dict["coreRadMag"]),
        minSegLen = convert(Float64, dict["minSegLen"]),
        maxSegLen = convert(Float64, dict["maxSegLen"]),
        minArea = convert(Float64, dict["minArea"]),
        maxArea = convert(Float64, dict["maxArea"]),
        maxConnect = convert(Int, dict["maxConnect"]),
        remesh = dict["remesh"],
        collision = dict["collision"],
        separation = dict["separation"],
        virtualRemesh = dict["virtualRemesh"],
        parCPU = dict["parCPU"],
        parGPU = dict["parGPU"],
        edgeDrag = convert(Float64, dict["edgeDrag"]),
        screwDrag = convert(Float64, dict["screwDrag"]),
        climbDrag = convert(Float64, dict["climbDrag"]),
        lineDrag = convert(Float64, dict["lineDrag"]),
        mobility = mobDict[dict["mobility"]],
    )

    return DislocationParams
end

"""
```
loadParams(
    fileDislocationParameters::AbstractString,
    fileMaterialParameters::AbstractString,
    fileIntegrationParameters::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)
```
Loads simulation parameters out of a dictionary loaded from a JSON file. Returns a tuple of variable types ([`DislocationParameters`](@ref), [`MaterialParameters`](@ref), [`IntegrationParameters`](@ref), [`SlipSystem`](@ref), [`DislocationLoop`](@ref)) or vectors of those types.
"""
function loadParams(
    fileDislocationParameters::AbstractString,
    fileMaterialParameters::AbstractString,
    fileIntegrationParameters::AbstractString,
    fileSlipSystem::AbstractString,
    fileDislocationLoop::AbstractString,
)
    # We use JSON arrays because it lets us dump a variable number of args into a single JSON file. To keep things gonsistent we use them always. Hence the indices here.
    dictDislocationParameters = load(fileDislocationParameters)
    DislocationParams = loadDislocationParameters(dictDislocationParameters)

    dictMaterialParameters = load(fileMaterialParameters)
    MaterialParams = loadMaterialParameters(dictMaterialParameters)

    dictIntegrationParameters = load(fileIntegrationParameters)
    IntegrationParams = loadIntegrationParameters(dictIntegrationParameters)

    dictSlipSystem = load(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem)

    # There can be multiple dislocations per simulation parameters.
    dictDislocationLoop = load(fileDislocationLoop)
    if typeof(dictDislocationLoop) <: AbstractArray
        dislocationLoop = zeros(DislocationLoop, length(dictDislocationLoop))
        for i in eachindex(dislocationLoop)
            dislocationLoop[i] = loadDislocationLoop(dictDislocationLoop[i], slipSystems)
        end
    else
        dislocationLoop = loadDislocationLoop(dictDislocationLoop, slipSystems)
    end

    return DislocationParams, MaterialParams, IntegrationParams, slipSystems, dislocationLoop
end

"""
```
loadNetwork(fileDislocationNetwork::AbstractString)
```
Loads a dislocation network from a JSON file. Returns a [`DislocationNetwork`](@ref).
"""
function loadNetwork(fileDislocationNetwork::AbstractString)
    dict = load(fileDislocationNetwork)

    lenLinks = length(dict["links"])
    lenCoord = length(dict["coord"])
    numNode = [convert(Int, dict["numNode"][1])]
    numSeg = [convert(Int, dict["numSeg"][1])]
    maxConnect = [convert(Int, dict["maxConnect"][1])]
    links = zeros(Int, 2, lenLinks)
    slipPlane = zeros(3, lenLinks)
    bVec = zeros(3, lenLinks)
    coord = zeros(3, lenCoord)
    connectivity = zeros(Int, 2 * maxConnect[1] + 1, lenLinks)
    linksConnect = zeros(Int, 2, lenLinks)
    segIdx = zeros(Int, lenLinks, 3)
    segForce = zeros(3, 2, lenLinks)
    nodeVel = zeros(3, lenLinks)
    nodeForce = zeros(3, lenLinks)

    for i in 1:lenLinks
        links[:, i] = dict["links"][i]
        linksConnect[:, i] = dict["linksConnect"][i]
        segForce[:, 1, i] = dict["segForce"][i][1]
        segForce[:, 2, i] = dict["segForce"][i][2]
    end
    for i in 1:lenCoord
        slipPlane[:, i] = dict["slipPlane"][i]
        bVec[:, i] = dict["bVec"][i]
        coord[:, i] = dict["coord"][i]
        nodeVel[:, i] = dict["nodeVel"][i]
        nodeForce[:, i] = dict["nodeForce"][i]
        connectivity[:, i] = dict["connectivity"][i]
    end

    for i in 1:3
        segIdx[:, i] = dict["segIdx"][i]
    end

    dislocationNetwork = DislocationNetwork(;
        links = links,
        slipPlane = slipPlane,
        bVec = bVec,
        coord = coord,
        label = nodeType.(dict["label"]),
        segForce = segForce,
        nodeVel = nodeVel,
        nodeForce = nodeForce,
        numNode = numNode,
        numSeg = numSeg,
        maxConnect = maxConnect,
        linksConnect = linksConnect,
        connectivity = connectivity,
        segIdx = segIdx,
    )

    return dislocationNetwork
end

function loadIntegrationTime(fileIntegrationTime::AbstractString)
    dict = load(fileIntegrationTime)
    integrationTime = IntegrationTime(;
        dt = convert(Float64, dict["dt"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int, dict["step"]),
    )
    return integrationTime
end

function loadIntegrationTime(dict::Dict{T1, T2}) where {T1, T2}
    integrationTime = IntegrationTime(;
        dt = convert(Float64, dict["dt"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int, dict["step"]),
    )
    return integrationTime
end
