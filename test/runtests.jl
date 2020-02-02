using SafeTestsets

@safetestset "IO" begin
    include("ioTests.jl")
end
@safetestset "Generate dln" begin
    include("genDlnTests.jl")
end
