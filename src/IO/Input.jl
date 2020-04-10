function load(filename::AbstractString)
    dict = JSON.parsefile(filename)
    return dict
end

function loadCSV(
    filename::AbstractString;
    header = 0,
    transpose::Bool = false,
    delim = ',',
)
    df = CSV.read(
        filename;
        copycols = false,
        header = header,
        transpose = transpose,
        delim = delim,
    )
    @inbounds for i in names(df)
        ismissing(df[1, i]) ? df[1, i] = 0 : nothing
    end
    return df
end

function loadSlipSys(dict::Dict{T1, T1}) where {T1, T2}

    crystalStruct = makeTypeDict(AbstractCrystalStruct)

    slipSystem = SlipSystem(
        name = crystalStruct[dict["crystalStruct"]],
        slipPlane["slipPlane"],
        bVec = dict["bVec"],
    )

    return slipSystem
end

function loadDln(df::DataFrame, slipSystems::AbstractArray{<:Real, N} where {N})
    nRow = nrow(df)
    sources = zeros(DislocationLoop, nRow)
    span = zeros(Float64, 2, 3)
    segTypes = makeTypeDict(AbstractDlnSeg)
    dlnTypes = makeTypeDict(AbstractDlnStr)
    dist = makeTypeDict(AbstractDistribution)

    @inbounds for i = 1:nRow
        st = split.(df[i, :segType], ";")
        segType = [segTypes[st[i]] for i = 1:length(st)]
        length(segType) == 1 ? segType = segType[1] : nothing
        sl = split.(df[i, :segLen], ";")
        segLen = parse.(Float64, sl)
        _slipSystem = df[i, :slipSystem]
        try
            ss = split.(df[i, :slipSystem], ";")
            _slipSystem = parse.(Int64, ss)
        catch err
        end
        lbl = split.(df[i, :label], ";")
        label = convert.(nodeType, parse.(Int64, lbl))
        spanmin = split.(df[i, :spanmin], ";")
        spanmax = split.(df[i, :spanmax], ";")
        span[1, :] .= parse.(Float64, spanmin)
        span[2, :] .= parse.(Float64, spanmax)
        sources[i] = DislocationLoop(
            loopType = dlnTypes[df[i, :loopType]],
            numSides = df[i, :numSides],
            nodeSide = convert(Int64, df[i, :nodeSide]),
            numLoops = convert(Int64, df[i, :numLoops]),
            segType = segType,
            segLen = segLen,
            slipSystem = _slipSystem,
            _slipPlane = slipSystems[_slipSystem, 1:3],
            _bVec = slipSystems[_slipSystem, 4:6],
            label = label,
            buffer = convert(Float64, df[i, :buffer]),
            range = span,
            dist = dist[df[1, :dist]],
        )
    end
    return sources
end

function dlnLoadParams(dict::Dict{T1, T2}) where {T1, T2}

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

    return dlnLoadParams
end

function matLoadParams(dict::Dict{T1, T2}) where {T1, T2}

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

function intLoadParams(dict::Dict{T1, T2}) where {T1, T2}

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
    dlnParams = dlnLoadParams(df)
    matParams = matLoadParams(df)
    intParams = intLoadParams(df)
    slipSystems = loadSlipSys(fileSlipSys)
    sources = loadDln(df2, slipSystems)
    return (dlnParams, matParams, intParams, slipSystems, sources)
end # function
