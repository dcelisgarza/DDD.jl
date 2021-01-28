"""
```
loadJSON(filename::AbstractString)
```
Wrapper for `JSON.parsefile(filename)`. Loads a `JSON` file as a dictionary.
"""
function loadJSON(filename::AbstractString)
    dict = JSON.parsefile(filename)
    return dict
end

"""
```
loadDislocationLoop(dict::Dict{T1,T2} where {T1,T2}, slipSystem::SlipSystem)
```
Constructs [`DislocationLoop`](@ref) out of a dictionary and [`SlipSystem`](@ref) structure.
"""
function loadDislocationLoop(dict::Dict{T1,T2} where {T1,T2}, slipSystem::SlipSystem)
    dlnTypes = makeTypeDict(AbstractDlnStr)
    distributions = makeTypeDict(AbstractDistribution)

    slipPlane = slipSystem.slipPlane
    bVec = slipSystem.bVec

    range = SMatrix{3,2}(convert.(Float64, vcat(dict["range"]...)))

    numSides = convert(Int, dict["numSides"])
    nodeSide = convert(Int, dict["nodeSide"])
    numLoops = convert(Int, dict["numLoops"])
    nodeLoop = numSides * nodeSide

    dislocationLoop = DislocationLoop(;
        loopType = dlnTypes[dict["loopType"]],
        numSides = numSides,
        nodeSide = nodeSide,
        numLoops = numLoops,
        segLen = SVector{length(dict["segLen"]),Float64}(dict["segLen"]),
        slipSystem = convert.(Int, dict["slipSystem"]),
        _slipPlane = convert.(Float64, slipPlane[:, dict["slipSystem"]]),
        _bVec = convert.(Float64, bVec[:, dict["slipSystem"]]),
        label = SVector{nodeLoop,nodeTypeDln}(vcat(dict["label"]...)),
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
Constructs [`MaterialParameters`](@ref) out of a dictionary.
"""
function loadMaterialParameters(dict::Dict{T1,T2}) where {T1,T2}
    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    MaterialParams = MaterialParameters(;
        crystalStruct = crystalStruct[dict["crystalStruct"]],
        μ = convert(Float64, dict["μ"]),
        μMag = convert(Float64, dict["μMag"]),
        ν = convert(Float64, dict["ν"]),
        σPN = convert(Float64, dict["σPN"]),
    )

    return MaterialParams
end

"""
```
loadFEMParameters(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`FEMParameters`](@ref) out of a dictionary.
"""
function loadFEMParameters(dict::Dict{T1,T2}) where {T1,T2}
    meshDict = makeTypeDict(AbstractMesh)
    orderDict = makeTypeDict(AbstractElementOrder)
    modelDict = makeTypeDict(AbstractModel)

    FemParams = FEMParameters(;
        type = meshDict[dict["type"]],
        order = orderDict[dict["order"]],
        model = modelDict[dict["model"]],
        dx = convert(Float64, dict["dx"]),
        dy = convert(Float64, dict["dy"]),
        dz = convert(Float64, dict["dz"]),
        mx = convert(Int, dict["mx"]),
        my = convert(Int, dict["my"]),
        mz = convert(Int, dict["mz"]),
    )

    return FemParams
end

"""
```
loadBoundaries(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`Boundaries`](@ref) out of a dictionary.
!!! note
    `tK` may be null if it was factorised when the variable was saved.
"""
function loadBoundaries(dict::Dict{T1,T2}) where {T1,T2}
    uGammaDict = dict["uGamma"]
    mGammaDict = dict["mGamma"]
    tGammaDict = dict["tGamma"]

    uGamma = try (type = nodeTypeFE.(uGammaDict["type"]),
                    idx = Int.(uGammaDict["idx"]),
                    node = Int.(uGammaDict["node"]))
    catch err
        uGamma = nothing
    end

    mGamma = try (type = nodeTypeFE.(mGammaDict["type"]),
                    idx = Int.(mGammaDict["idx"]),
                    node = Int.(mGammaDict["node"]))
    catch err
        mGamma = nothing
    end

    tGamma = try (type = nodeTypeFE.(tGammaDict["type"]),
                    idx = Int.(tGammaDict["idx"]),
                    node = Int.(tGammaDict["node"]))
    catch err
        tGamma = nothing
    end

    uDofs = try Int.(dict["uDofs"]); catch err; nothing end
    mDofs = try Int.(dict["mDofs"]); catch err; nothing end
    tDofs = try Int.(dict["tDofs"]); catch err; nothing end

    typeof(dict["tK"]) <: AbstractArray ? tK = Float64.(dict["tK"]) : tK = nothing

    return Boundaries(; 
            uGamma = uGamma,
            tGamma = tGamma,
            mGamma = mGamma,
            uDofs = uDofs,
            tDofs = tDofs,
            mDofs = mDofs,
            tK = tK
        )
end

"""
```
loadForceDisplacement(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`ForceDisplacement`](@ref) out of a dictionary. It makes the arrays sparse and drops zeros under `eps(Float64)`.
"""
function loadForceDisplacement(dict::Dict{T1,T2}) where {T1,T2}
    return ForceDisplacement(;
        uTilde = droptol!(sparse(Float64.(dict["uTilde"])), eps()),
        uHat = droptol!(sparse(Float64.(dict["uHat"])), eps()),
        u = droptol!(sparse(Float64.(dict["u"])), eps()),
        fTilde = droptol!(sparse(Float64.(dict["fTilde"])), eps()),
        fHat = droptol!(sparse(Float64.(dict["fHat"])), eps()),
        f = droptol!(sparse(Float64.(dict["f"])), eps()),
    )
end

"""
```
loadIntegrationParameters(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`IntegrationParameters`](@ref) out of a dictionary.
"""
function loadIntegrationParameters(dict::Dict{T1,T2}) where {T1,T2}
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
loadSlipSystem(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`SlipSystem`](@ref) out of a dictionary.
"""
function loadSlipSystem(dict::Dict{T1,T2}) where {T1,T2}
    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    lenSlipSys = length(dict["slipPlane"])

    slipPlane = SMatrix{3,lenSlipSys}(convert.(Float64, vcat(dict["slipPlane"]...)))
    bVec = SMatrix{3,lenSlipSys}(convert.(Float64, vcat(dict["bVec"]...)))

    slipSystem = SlipSystem(;
        crystalStruct = crystalStruct[dict["crystalStruct"]],
        slipPlane = slipPlane,
        bVec = bVec,
    )

    return slipSystem
end

"""
```
loadDislocationParameters(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`DislocationParameters`](@ref) out of a dictionary.
"""
function loadDislocationParameters(dict::Dict{T1,T2}) where {T1,T2}
    mobDict = makeTypeDict(AbstractMobility)

    dragCoeffs = dict["dragCoeffs"]
    
    DislocationParams = DislocationParameters(;
        coreRad = convert(Float64, dict["coreRad"]),
        dragCoeffs = namedtuple(dict["dragCoeffs"]),
        coreRadMag = convert(Float64, dict["coreRadMag"]),
        coreEnergy = convert(Float64, dict["coreEnergy"]),
        minSegLen = convert(Float64, dict["minSegLen"]),
        maxSegLen = convert(Float64, dict["maxSegLen"]),
        minArea = convert(Float64, dict["minArea"]),
        maxArea = convert(Float64, dict["maxArea"]),
        remesh = dict["remesh"],
        collision = dict["collision"],
        separation = dict["separation"],
        virtualRemesh = dict["virtualRemesh"],
        parCPU = dict["parCPU"],
        parGPU = dict["parGPU"],
        mobility = mobDict[dict["mobility"]],
    )

    return DislocationParams
end

"""
```
loadNetwork(fileDislocationNetwork::AbstractString)
```
Constructs [`DislocationNetwork`](@ref) from a `JSON` file.
"""
function loadNetwork(fileDislocationNetwork::AbstractString)
    dict = loadJSON(fileDislocationNetwork)

    lenLinks = length(dict["links"])
    lenCoord = length(dict["coord"])
    numNode = [convert(Int, dict["numNode"][1])]
    numSeg = [convert(Int, dict["numSeg"][1])]
    maxConnect = convert(Int, dict["maxConnect"])
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
        label = nodeTypeDln.(dict["label"]),
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
"""
```
loadNetwork(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`DislocationNetwork`](@ref) from a dictionary.
"""
function loadNetwork(dict::Dict{T1,T2}) where {T1,T2}
    lenLinks = length(dict["links"])
    lenCoord = length(dict["coord"])
    numNode = [convert(Int, dict["numNode"][1])]
    numSeg = [convert(Int, dict["numSeg"][1])]
    maxConnect = convert(Int, dict["maxConnect"])
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
        label = nodeTypeDln.(dict["label"]),
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

"""
```
loadIntegrationTime(fileIntegrationTime::AbstractString)
```
Constructs [`IntegrationTime`](@ref) from a `JSON` file.
"""
function loadIntegrationTime(fileIntegrationTime::AbstractString)
    dict = loadJSON(fileIntegrationTime)
    integrationTime = IntegrationTime(;
        dt = convert(Float64, dict["dt"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int, dict["step"]),
    )
    return integrationTime
end
"""
```
loadIntegrationTime(dict::Dict{T1,T2}) where {T1,T2}
```
Constructs [`IntegrationTime`](@ref) from a dictionary.
"""
function loadIntegrationTime(dict::Dict{T1,T2}) where {T1,T2}
    integrationTime = IntegrationTime(;
        dt = convert(Float64, dict["dt"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int, dict["step"]),
    )
    return integrationTime
end

"""
```
loadParameters(
    fileDislocationParameters::T,
    fileMaterialParameters::T,
    fileFEMParameters::T,
    fileIntegrationParameters::T,
    fileSlipSystem::T,
    fileDislocationLoop::T,
) where {T <: AbstractString}
```
Constructs simulation parameters, ([`DislocationParameters`](@ref), [`MaterialParameters`](@ref), [`FEMParameters`](@ref), [`IntegrationParameters`](@ref), [`SlipSystem`](@ref), [`DislocationLoop`](@ref)) from `JSON` files.
"""
function loadParameters(
    fileDislocationParameters::T,
    fileMaterialParameters::T,
    fileFEMParameters::T,
    fileIntegrationParameters::T,
    fileSlipSystem::T,
    fileDislocationLoop::T,
) where {T <: AbstractString}
    # We use JSON arrays because it lets us dump a variable number of args into a single JSON file. To keep things gonsistent we use them always. Hence the indices here.
    dictDislocationParameters = loadJSON(fileDislocationParameters)
    DislocationParams = loadDislocationParameters(dictDislocationParameters)

    dictMaterialParameters = loadJSON(fileMaterialParameters)
    MaterialParams = loadMaterialParameters(dictMaterialParameters)

    dictFEMParameters = loadJSON(fileFEMParameters)
    FemParams = loadFEMParameters(dictFEMParameters)

    dictIntegrationParameters = loadJSON(fileIntegrationParameters)
    IntegrationParams = loadIntegrationParameters(dictIntegrationParameters)

    dictSlipSystem = loadJSON(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem)

    # There can be multiple dislocations per simulation parameters.
    dictDislocationLoop = loadJSON(fileDislocationLoop)
    if typeof(dictDislocationLoop) <: AbstractArray
        dislocationLoop = zeros(DislocationLoop, length(dictDislocationLoop))
        for i in eachindex(dislocationLoop)
            dislocationLoop[i] =
                loadDislocationLoop(dictDislocationLoop[i], slipSystems)
        end
    else
        dislocationLoop = loadDislocationLoop(dictDislocationLoop, slipSystems)
    end

    return DislocationParams,
    MaterialParams,
    FemParams,
    IntegrationParams,
    slipSystems,
    dislocationLoop
end