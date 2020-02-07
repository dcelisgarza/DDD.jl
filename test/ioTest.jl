using DDD
using Test, DataFrames

cd(@__DIR__)
# Filenames
inFilename = "../inputs/simParams/sampleParams.csv"
outFilename = "../outputs/simParams/sampleParams.csv"

# Load parameters.
dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)

# Testing inputs
@testset "Compare loading inputs and outputs" begin
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end

@testset "Clean DataFrame" begin
    df = DataFrame(numSources = Any[], fieldName = Any[])
    push!(df, ([1, 2], [2, 4]))
    @test cleanFieldDf(df, :fieldName, Integer) == [2; 4]

    df = DataFrame(numSources = Any[], fieldName = Any[])
    push!(df, (1, 2))
    push!(df, (2, 4))
    @test cleanFieldDf(df, :fieldName, Integer) == [2; 4]
end
