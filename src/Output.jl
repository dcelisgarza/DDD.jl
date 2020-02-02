module Output
using CSV, DataFrames
# using ..CustomTypes
using ..DislocationBase
using ..MaterialBase
using ..CustomIntegration
using ..DdFemBase

export pushToDataFrame!, saveParams

"""
Pushes data to a dataframe for saving later.
"""
function pushToDataFrame!(
    df::DataFrame,
    data::Union{DislocationP,MaterialP,IntegrationP},
)
    fieldNames = fieldnames(typeof(data))
    # lenFieldNames = length(fieldNames)
    for i in eachindex(fieldNames)
        push!(df, (fieldNames[i], getproperty(data, fieldNames[i])))
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
    filename::AbstractString = "out",
    extension::AbstractString = ".csv";
    delim::Char = ',',
)
    df = DataFrame(var = Any[], val = Any[])

    pushToDataFrame!(df, dlnParams)
    pushToDataFrame!(df, matParams)
    pushToDataFrame!(df, intParams)

    CSV.write(filename*"Params"*extension, df; delim = delim, writeheader = false)
end

end # Output
