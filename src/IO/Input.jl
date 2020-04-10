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

function loadSlipSys(filename::AbstractString, delim = ',')
    return readdlm(filename, delim)
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

function dlnLoadParams(df::DataFrame)
    mobDict = makeTypeDict(AbstractMobility)
    dlnLoadParams = DislocationP(
        coreRad = convert(Float64, df[1, :coreRad]),
        coreRadMag = convert(Float64, df[1, :coreRadMag]),
        minSegLen = convert(Float64, df[1, :minSegLen]),
        maxSegLen = convert(Float64, df[1, :maxSegLen]),
        minArea = convert(Float64, df[1, :minArea]),
        maxArea = convert(Float64, df[1, :maxArea]),
        maxConnect = convert(Int64, df[1, :maxConnect]),
        remesh = df[1, :remesh],
        collision = df[1, :collision],
        separation = df[1, :separation],
        virtualRemesh = df[1, :virtualRemesh],
        edgeDrag = convert(Float64, df[1, :edgeDrag]),
        screwDrag = convert(Float64, df[1, :screwDrag]),
        climbDrag = convert(Float64, df[1, :climbDrag]),
        lineDrag = convert(Float64, df[1, :lineDrag]),
        mobility = mobDict[df[1, :mobility]],
    )
    return dlnLoadParams
end

function matLoadParams(df::DataFrame)
    strucDict = makeTypeDict(AbstractCrystalStruct)
    matParams = MaterialP(
        μ = convert(Float64, df[1, :μ]),
        μMag = convert(Float64, df[1, :μMag]),
        ν = convert(Float64, df[1, :ν]),
        E = convert(Float64, df[1, :E]),
        crystalStruct = strucDict[df[1, :crystalStruct]],
    )
    return matParams
end

function intLoadParams(df::DataFrame)
    integDict = makeTypeDict(AbstractIntegrator)
    intParams = IntegrationP(
        dt = convert(Float64, df[1, :dt]),
        tmin = convert(Float64, df[1, :tmin]),
        tmax = convert(Float64, df[1, :tmax]),
        method = integDict[df[1, :method]],
        abstol = convert(Float64, df[1, :abstol]),
        reltol = convert(Float64, df[1, :reltol]),
        time = convert(Float64, df[1, :time]),
        step = convert(Int64, df[1, :step]),
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
