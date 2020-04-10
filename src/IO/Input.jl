function load(filename::AbstractString)
    dict = JSON.parsefile(filename)
    return dict
end

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

function loadDislocationLoop(
    dict::Dict{T1, T2} where {T1, T2},
    slipSystem::SlipSystem,
)
    dlnTypes = makeTypeDict(AbstractDlnStr)
    distributions = makeTypeDict(AbstractDistribution)

    slipPlane = slipSystem.slipPlane
    bVec = slipSystem.bVec

    dislocationLoop = DislocationLoop(
        loopType = dlnTypes[dict["loopType"]],
        numSides = convert(Int64, dict["numSides"]),
        nodeSide = convert(Int64, dict["nodeSide"]),
        numLoops = convert(Int64, dict["numLoops"]),
        segLen = convert.(Float64, dict["segLen"]),
        slipSystem = convert(Int64, dict["slipSystem"]),
        _slipPlane = convert.(Float64, slipPlane[dict["slipSystem"], 1:3]),
        _bVec = convert.(Float64, bVec[dict["slipSystem"], 1:3]),
        label = nodeType.(dict["label"]),
        buffer = convert.(Float64, dict["buffer"]),
        range = [0.0 0.0 0.0; 0.0 0.0 0.0],
        dist = distributions[dict["dist"]],
    )

    return dislocationLoop
end

function loadDislocationP(dict::Dict{T1, T2}) where {T1, T2}

    mobDict = makeTypeDict(AbstractMobility)

    dlnParams = DislocationP(
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

    return dlnParams
end

function loadMaterialP(dict::Dict{T1, T2}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    matParams = MaterialP(
        μ = convert(Float64, dict["μ"]),
        μMag = convert(Float64, dict["μMag"]),
        ν = convert(Float64, dict["ν"]),
        E = convert(Float64, dict["E"]),
        crystalStruct = crystalStruct[dict["crystalStruct"]],
    )

    return matParams
end

function loadIntegrationP(dict::Dict{T1, T2}) where {T1, T2}

    integDict = makeTypeDict(AbstractIntegrator)

    intParams = IntegrationP(
        dt = convert(Float64, dict["dt"]),
        tmin = convert(Float64, dict["tmin"]),
        tmax = convert(Float64, dict["tmax"]),
        method = integDict[dict["method"]],
        abstol = convert(Float64, dict["abstol"]),
        reltol = convert(Float64, dict["reltol"]),
        time = convert(Float64, dict["time"]),
        step = convert(Int64, dict["step"]),
    )

    return intParams
end

function loadParams(
    fileParams::AbstractString,
    fileSlipSys::AbstractString,
    fileDln::AbstractString,
)
    df = loadCSV(
        fileParams::AbstractString;
        header = 1,
        transpose = true,
        delim = ',',
    )
    df2 = loadCSV(
        fileDln::AbstractString;
        header = 1,
        transpose = true,
        delim = ',',
    )
    dlnParams = loadDislocationP(df)
    matParams = loadMaterialP(df)
    intParams = loadIntegrationP(df)
    slipSystems = loadSlipSystem(fileSlipSys)
    sources = loadDislocationLoop(df2, slipSystems)
    return (dlnParams, matParams, intParams, slipSystems, sources)
end # function
