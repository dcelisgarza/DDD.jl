JSON.lower(t::T) where {T<:Union{AbstractCrystalStruct, AbstractMobility, AbstractIntegrator, AbstractDlnSeg, AbstractDlnStr, AbstractDistribution}} = string(t)
JSON.lower(t::nodeType) = Int(t)

function save(filename::AbstractString, args...)
    open(filename, "w") do io
        JSON.print(io, args)
    end
end

"""
Pushes data to a dataframe for saving later.
"""
function pushToDataFrame!(
    df::DataFrame,
    data::Union{DislocationP, MaterialP, IntegrationP},
)
    fieldNames = fieldnames(typeof(data))
    for fieldName in fieldNames
        push!(df, (fieldName, getproperty(data, fieldName)))
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
