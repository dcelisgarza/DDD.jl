"""
Pushes data to a dataframe for saving later.
"""
function pushToDataFrame!(
    df::DataFrame,
    data::Union{DislocationP, MaterialP, IntegrationP},
)
    fieldNames = fieldnames(typeof(data))
    for fieldName in fieldNames
        push!(df, (fieldNames, getproperty(data, fieldNames)))
    end
    return df
end # pushToDataFrame
"""
Saves simulation parameters.
"""
function saveParams(
    dlnParams::DislocationP,
    matParams::MaterialP,
    intParams::IntegrationP,
    filename::AbstractString;
    delim::Char = ',',
)
    df = DataFrame(var = Any[], val = Any[])

    pushToDataFrame!(df, dlnParams)
    pushToDataFrame!(df, matParams)
    pushToDataFrame!(df, intParams)

    CSV.write(filename, df; delim = delim, writeheader = false)
end
