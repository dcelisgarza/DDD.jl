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
    for i in names(df)
        ismissing(df[1, i]) ? df[1, i] = 0 : nothing
    end
    return df
end

function loadSlipSys(filename::AbstractString, delim = ',')
    slipSystems = readdlm(filename, delim)
    return slipSystems
end

function loadDln(df::DataFrame, slipSystems::AbstractArray{<:Real,N} where {N})
    nRow = nrow(df)
    sources = zeros(DislocationLoop, nRow)
    span = zeros(2, 3)
    segTypes = Dict(
        "segEdge()" => segEdge(),
        "segEdgeN()" => segEdgeN(),
        "segScrew()" => segScrew(),
        "segMixed()" => segMixed(),
        "segNone()" => segNone(),
        "DDD.segEdge()" => segEdge(),
        "DDD.segEdgeN()" => segEdgeN(),
        "DDD.segScrew()" => segScrew(),
        "DDD.segMixed()" => segMixed(),
        "DDD.segNone()" => segNone(),
    )
    dist = Dict(
        "Zeros()" => Zeros(),
        "Rand()" => Rand(),
        "Randn()" => Randn(),
        "Regular()" => Regular(),
        "DDD.Zeros()" => Zeros(),
        "DDD.Rand()" => Rand(),
        "DDD.Randn()" => Randn(),
        "DDD.Regular()" => Regular(),
    )
    for i = 1:nRow
        st = split.(df[i, :segType], ";")
        segType = [segTypes[st[i]] for i = 1:length(st)]
        sl = split.(df[i, :segLen], ";")
        segLen = parse.(Float64, sl)
        ss = split.(df[i, :slipSystem], ";")
        slipSystem = parse.(Int64, ss)
        lbl = split.(df[i, :label], ";")
        label = convert.(nodeType, parse.(Int64, lbl))
        spanmin = split.(df[i, :spanmin], ";")
        spanmax = split.(df[i, :spanmax], ";")
        span[1, :] .= parse.(Float64, spanmin)
        span[2, :] .= parse.(Float64, spanmax)
        sources[i] = DislocationLoop(
            loopSides(df[i, :numSides]),
            df[i, :nodeSide],
            df[i, :numLoops],
            segType,
            segLen,
            slipSystem,
            slipSystems[slipSystem, 1:3],
            slipSystems[slipSystem, 4:6],
            label,
            convert(Float64,df[i, :buffer]),
            span,
            dist[df[1, :dist]],
        )
    end
    return sources
end

function dlnLoadParams(df::DataFrame)
    mobDict = Dict(
        "mobBCC()" => mobBCC(),
        "mobFCC()" => mobFCC(),
        "mobHCP()" => mobHCP(),
        "DDD.mobBCC()" => mobBCC(),
        "DDD.mobFCC()" => mobFCC(),
        "DDD.mobHCP()" => mobHCP(),
    )
    dlnLoadParams = DislocationP(
        df[1, :coreRad],
        df[1, :coreRadMag],
        df[1, :minSegLen],
        df[1, :maxSegLen],
        df[1, :minArea],
        df[1, :maxArea],
        convert(Int64, df[1, :maxConnect]),
        df[1, :remesh],
        df[1, :collision],
        df[1, :separation],
        df[1, :virtualRemesh],
        df[1, :edgeDrag],
        df[1, :screwDrag],
        df[1, :climbDrag],
        df[1, :lineDrag],
        mobDict[df[1, :mobility]],
    )
    return dlnLoadParams
end

function matLoadParams(df::DataFrame)
    strucDict = Dict(
        "BCC()" => BCC(),
        "FCC()" => FCC(),
        "HCP()" => HCP(),
        "DDD.BCC()" => BCC(),
        "DDD.FCC()" => FCC(),
        "DDD.HCP()" => HCP(),
    )
    matParams = MaterialP(
        df[1, :μ],
        df[1, :μMag],
        df[1, :ν],
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
        df[1, :dt],
        df[1, :tmin],
        df[1, :tmax],
        integDict[df[1, :method]],
        df[1, :abstol],
        df[1, :reltol],
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
