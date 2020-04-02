using SafeTestsets

@safetestset "IO" begin
    include("./ioTest.jl")
end
@safetestset "Generate dln" begin
    include("./dlnNetworkTest.jl")
end
@safetestset "FEM" begin
    include("./FEMTest.jl")
end
@safetestset "Misc tests" begin
    include("./MiscTest.jl")
end
@safetestset "Plotting" begin
    include("PlotTest.jl")
end
# @safetestset "Performant code" begin
#     include("performantCodeTest.jl")
# end
