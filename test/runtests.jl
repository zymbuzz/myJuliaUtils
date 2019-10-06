println("Testing...")

using myJuliaUtils
using Test

@time @testset "Testing vech" begin include("TestingVech.jl") end
@time @testset "Testing Wish" begin include("TestWish/TestingWish.jl") end
@time @testset "Testing cholx and randx" begin include("TestingCholXRandX.jl") end
@time @testset "Testing some functions in extfunc" begin include("TestingExtfunc.jl") end
@time @testset "Testing ACF function" begin include("testAcf/testAcf.jl") end
@time @testset "Testing decVCV function" begin include("testdecVCV.jl") end
@time @testset "Testing the fitting of distributions" begin include("testFitDist.jl") end
