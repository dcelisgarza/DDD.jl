using DDD
using Test

import DDD: compStruct

inFilename = "../inputs/simParams/sample"
outFilename = "../outputs/simParams/sample"
dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)


@testset "DDD.jl" begin
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end