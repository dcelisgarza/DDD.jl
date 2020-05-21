using SafeTestsets

@safetestset "Miscelaneous" begin
    include("./MiscTest.jl")
end
@safetestset "Generate dln" begin
    include("./dlnNetworkTest.jl")
end
@safetestset "IO" begin
    include("./ioTest.jl")
end
@safetestset "Topology" begin
    include("./TopologyTest.jl")
end
@safetestset "FEM" begin
    include("./FEMTest.jl")
end
@safetestset "Post-process" begin
    include("./PlotTest.jl")
end
# @safetestset "Performant code" begin
#     include("performantCodeTest.jl")
# end
