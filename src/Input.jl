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
            dlnTypes[df[i, :loopType]],
            df[i, :numSides],
            convert(Int64, df[i, :nodeSide]),
            convert(Int64, df[i, :numLoops]),
            segType,
            segLen,
            _slipSystem,
            slipSystems[_slipSystem, 1:3],
            slipSystems[_slipSystem, 4:6],
            label,
            convert(Float64, df[i, :buffer]),
            span,
            dist[df[1, :dist]],
        )
    end
    return sources
end

function dlnLoadParams(df::DataFrame)
    mobDict = makeTypeDict(AbstractMobility)
    dlnLoadParams = DislocationP(
        convert(Float64, df[1, :coreRad]),
        convert(Float64, df[1, :coreRadMag]),
        convert(Float64, df[1, :minSegLen]),
        convert(Float64, df[1, :maxSegLen]),
        convert(Float64, df[1, :minArea]),
        convert(Float64, df[1, :maxArea]),
        convert(Int64, df[1, :maxConnect]),
        df[1, :remesh],
        df[1, :collision],
        df[1, :separation],
        df[1, :virtualRemesh],
        convert(Float64, df[1, :edgeDrag]),
        convert(Float64, df[1, :screwDrag]),
        convert(Float64, df[1, :climbDrag]),
        convert(Float64, df[1, :lineDrag]),
        mobDict[df[1, :mobility]],
    )
    return dlnLoadParams
end

function matLoadParams(df::DataFrame)
    strucDict = makeTypeDict(AbstractCrystalStruct)
    matParams = MaterialP(
        convert(Float64, df[1, :μ]),
        convert(Float64, df[1, :μMag]),
        convert(Float64, df[1, :ν]),
        strucDict[df[1, :crystalStruct]],
    )
    return matParams
end

function intLoadParams(df::DataFrame)
    integDict = Dict(
        "CustomTrapezoid()" => CustomTrapezoid(),
        "DDD.CustomTrapezoid()" => CustomTrapezoid(),
    )
    intParams = IntegrationP(
        convert(Float64, df[1, :dt]),
        convert(Float64, df[1, :tmin]),
        convert(Float64, df[1, :tmax]),
        integDict[df[1, :method]],
        convert(Float64, df[1, :abstol]),
        convert(Float64, df[1, :reltol]),
        convert(Float64, df[1, :time]),
        convert(Int64, df[1, :step]),
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
    return dlnParams, matParams, intParams, slipSystems, sources
end # function
