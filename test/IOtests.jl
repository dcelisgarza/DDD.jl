using DDD
using Test

import DDD: compStruct

cd(@__DIR__)
inFilename = "../inputs/simParams/sampleParams.csv"
outFilename = "../outputs/simParams/sampleParams.csv"

dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)

@testset "Compare Parameter Structures" begin
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end
