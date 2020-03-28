using DDD
using Test, DataFrames

cd(@__DIR__)
@testset "Compare loading input and output parameters" begin
    # Filenames
    inFileParams = "../inputs/simParams/sampleParams.csv"
    inFileSlipSys = "../data/slipSystems/bcc.csv"
    inFileDln = "../inputs/dln/samplePrismaticShear.csv"
    outFileParams = "../outputs/simParams/sampleParams.csv"
    # Load parameters.
    dlnParams, matParams, intParams, slipSystems, sources =
        loadParams(inFileParams, inFileSlipSys, inFileDln)
    saveParams(dlnParams, matParams, intParams, outFileParams; delim = ',')
    dlnParams2, matParams2, intParams2, slipSystems, sources =
        loadParams(outFileParams, inFileSlipSys, inFileDln)
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end
