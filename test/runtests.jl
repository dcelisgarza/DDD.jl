using SafeTestsets

@safetestset "Miscelaneous" begin
    include("./MiscTest.jl")
end
@safetestset "IO" begin
    include("./ioTest.jl")
end
@safetestset "Construct Dln" begin
    include("./ConstructNetworkTest.jl")
end
@safetestset "Segment Forces" begin
    include("./SegmentForcesTest.jl")
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
