using DDD
using Test

inFilename = "C:/Users/Daniel Celis Garza/.julia/dev/DDD/inputs/simParams/sample"
outFilename = "C:/Users/Daniel Celis Garza/.julia/dev/DDD/outputs/simParams/sample"
dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
dlnParams2, matParams2, intParams2 = loadParams(outFilename)


function compare(arg1, arg2)
    @assert typeof(arg1) == typeof(arg2)

    names = fieldnames(typeof(arg1))
    for i in names
        result = getproperty(arg1, i) == getproperty(arg1, i)
        if result == false return false end
    end
    return true
end

@testset "DDD.jl" begin
    @test compare(dlnParams, dlnParams2)
    @test compare(matParams, matParams2)
    @test compare(intParams, intParams2)
end
