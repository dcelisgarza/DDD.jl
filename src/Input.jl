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

function loadInitDln(df::DataFrame, slipSystems::AbstractArray{<:Real,N} where {N})
    nRow = nrow(df)
    sources = zeros(DislocationLoop, difLoops)
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
    for i = 1:nRow
        st = split.(df[i, :segType], ";")
        segType = [dict[st[i]] for i = 1:length(st)]
        sl = split.(df[i, :segLen], ";")
        segLen = parse.(Float64, sl)
        ss = split.(df[i, :slipSystem], ";")
        slipSystem = parse.(Int, ss)
        lbl = split.(df[i,:label],";")
        label = convert.(nodeType,parse.(Int, lbl))
        sources[i] = DislocationLoop(
            loopSides(df[i, :numSides]),
            df[i, :nodeSide],
            segType,
            segLen,
            slipSystems[slipSystem, 1:3],
            slipSystems[slipSystem, 4:6],
            label,
            df[i,:numLoops]
        )
    end
    return sources
end

function dlnLoadParams(df::DataFrame)
    mobDict = Dict("mobBCC()" => mobBCC(), "mobFCC()" => mobFCC(), "mobHCP()" => mobHCP(),
    "DDD.mobBCC()" => mobBCC(), "DDD.mobFCC()" => mobFCC(), "DDD.mobHCP()" => mobHCP())
    dlnLoadParams = DislocationP(
        df[1, :coreRad],
        df[1, :coreRadMag],
        df[1, :minSegLen],
        df[1, :maxSegLen],
        df[1, :minArea],
        df[1, :maxArea],
        convert(Integer, df[1, :maxConnect]),
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
    strucDict = Dict("BCC()" => BCC(), "FCC()" => FCC(), "HCP()" => HCP(),
    "DDD.BCC()" => BCC(), "DDD.FCC()" => FCC(), "DDD.HCP()" => HCP())
    matParams = MaterialP(
        df[1, :μ],
        df[1, :μMag],
        df[1, :ν],
        strucDict[df[1, :crystalStruct]],
    )
    return matParams
end

function intLoadParams(df::DataFrame)
    integDict = Dict("CustomTrapezoid()" => CustomTrapezoid(),
    "DDD.CustomTrapezoid()" => CustomTrapezoid())
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

function loadParams(filename::AbstractString,)
    df = loadCSV(
        filename::AbstractString;
        header = 1,
        transpose = true,
        delim = ',',
    )
    dlnParams = dlnLoadParams(df)
    matParams = matLoadParams(df)
    intParams = intLoadParams(df)
    return dlnParams, matParams, intParams
end # function
