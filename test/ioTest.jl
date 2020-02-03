using DDD
using Test

cd(@__DIR__)
# Filenames
inFilename = "../inputs/simParams/sampleParams.csv"
outFilename = "../outputs/simParams/sampleParams.csv"

# Load parameters.
dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)

# Testing inputs
@testset "Compare Parameter Structures" begin
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end
