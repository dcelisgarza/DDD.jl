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

# function cleanFieldDf(df::DataFrame, fieldname::Symbol, type::Type)
#     result = 0
#     if length(df[!, :numSources]) > 1
#         try
#             result = df[!, fieldname]
#         catch e
#             result = eval(Meta.parse(df[1, fieldname]))
#         end
#     else
#         result = df[1, fieldname]
#     end
#     return convert.(type, result)
# end

function dlnLoadParams(df::DataFrame)
    # numSources = cleanFieldDf(df, :numSources, Integer)
    # slipSystems = cleanFieldDf(df, :slipSystems, Integer)
    # distSource = cleanFieldDf(df, :distSource, Float64)

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
    )
    return dlnLoadParams
end

function matLoadParams(df::DataFrame)
    matParams = MaterialP(
        df[1, :μ],
        df[1, :μMag],
        df[1, :ν],
        df[1, :crystalStruct],
    )
    return matParams
end

function intLoadParams(df::DataFrame)
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
