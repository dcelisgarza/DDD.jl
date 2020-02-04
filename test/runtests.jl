using SafeTestsets

@safetestset "IO" begin
    include("./ioTest.jl")
end
@safetestset "Generate dln" begin
    include("./dlnNetworkTest.jl")
end
# @safetestset "Performant code" begin
#     include("performantCodeTest.jl")
# end
