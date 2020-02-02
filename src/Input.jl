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

function dlnLoadParams(df)
    numSources = 0
    if length(df[!, :numSources]) > 1
        try
            numSources = convert.(Integer, df[!, :numSources])
        catch e
            numSources = eval(Meta.parse(df[1, :numSources]))
        end
    else
        numSources = convert(Integer, df[1, :numSources])
    end

    slipSystems = 0
    if length(df[!, :slipSystems]) > 1
        try
            slipSystems = convert.(Integer, df[!, :slipSystems])
        catch e
            slipSystems = eval(Meta.parse(df[1, :slipSystems]))
        end
    else
        slipSystems = convert(Integer, df[1, :slipSystems])
    end

    distSource = 0
    if length(df[!, :slipSystems]) > 1
        try
            distSource = convert.(Float64, df[!, :distSource])
        catch e
            distSource = eval(Meta.parse(df[1, :distSource]))
        end
    else
        distSource = convert(Float64, df[1, :distSource])
    end

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
        df[1, :mobility],
        numSources,
        slipSystems,
        distSource,
        # convert.(Integer, df[!, :numSources]),
        # convert.(Integer, df[!, :slipSystems]),
        # convert.(Float64, df[!, :distSource]),
    )
    return dlnLoadParams
end

function matLoadParams(df)
    matParams = MaterialP(
        df[1, :μ],
        df[1, :μMag],
        df[1, :ν],
        df[1, :crystalStruct],
    )
    return matParams
end

function intLoadParams(df)
    intParams = IntegrationP(
        df[1, :dt],
        df[1, :tmin],
        df[1, :tmax],
        df[1, :method],
        df[1, :abstol],
        df[1, :reltol],
    )
    return intParams
end

function loadParams(
    filename::AbstractString,
    extension::AbstractString = ".csv",
)
    df = loadCSV(
        filename * "Params" * extension;
        header = 1,
        transpose = true,
        delim = ',',
    )
    dlnParams = dlnLoadParams(df)
    matParams = matLoadParams(df)
    intParams = intLoadParams(df)
    return dlnParams, matParams, intParams
end # function

function loadParams(
    filename::AbstractString,
    ::DislocationP,
    extension::AbstractString = ".csv",
)
    df = loadCSV(filename * "DlnP" * extension; header = 1, transpose = true)
    dlnParams = dlnLoadParams(df)
    return dlnParams
end # function

function loadParams(
    filename::AbstractString,
    ::MaterialP,
    extension::AbstractString = ".csv",
)
    df = loadCSV(filename * "MatP" * extension; header = 1, transpose = true)
    matParams = matLoadParams(df)
    return matParams
end # function

function loadParams(
    filename::AbstractString,
    ::IntegrationP,
    extension::AbstractString = ".csv",
)
    df = loadCSV(filename * "IntP" * extension; header = 1, transpose = true)
    intParams = IntegrationP(df)
    return intParams
end # function

function genSourcesBCC(
    dlnParams::DislocationP,
    filename = "D:/Projects/DDD/data/slipSystems/bcc.csv",
)
    df = loadCSV(filename; header = 1, transpose = false)
    dlnNetwork = zero(DislocationNetwork)


    return dlnNetwork

end # function
