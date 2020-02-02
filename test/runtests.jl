using SafeTestsets

@safetestset "IO" begin
    include("./ioTest.jl")
end
@safetestset "Generate dln" begin
    include("./genSegTest.jl")
end
