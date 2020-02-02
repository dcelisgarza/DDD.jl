using DDD
using Test

inFilename = "../inputs/simParams/sample"
outFilename = "../outputs/simParams/sample"
dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)

import DDD: compStruct
@testset "DDD.jl" begin
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end


# C:\Users\Daniel Celis Garza\.julia\dev\DDD\test\runtests.jl
# C:\Users\Daniel Celis Garza\.julia\dev\DDD\outputs\simParams\sampleParams.csv
