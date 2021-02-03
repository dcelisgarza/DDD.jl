using SafeTestsets

# 
# @time @safetestset "Miscelaneous" begin
#     include("./MiscTest.jl")
# end

# ##
# @time @safetestset "IO" begin
#     include("./ioTest.jl")
# end

# ##
# @time @safetestset "Construct Dln" begin
#     include("./ConstructNetworkTest.jl")
# end

# ##
# @time @safetestset "Segment Force" begin
#     include("./SegmentForcesTest.jl")
# end

# ##
# @time @safetestset "Mobility" begin
#     include("./MobilityTest.jl")
# end

# ##
# @time @safetestset "Topology" begin
#     include("./TopologyTest.jl")
# end

##
@time @safetestset "FEM" begin
    include("./FEMTest.jl")
end

##
@time @safetestset "Dislocation-FEM" begin
    include("./DislocationFEMTest.jl")
end

##
@time @safetestset "Post-process" begin
    include("./PlotTest.jl")
end

##
@time @safetestset "Integral tests" begin
    include("./IntegralTest.jl")
end