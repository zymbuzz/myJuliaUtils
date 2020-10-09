println("Testing...")

using myJuliaUtils
using Test
using SafeTestsets

@time @safetestset "Testing vech" begin include("TestingVech.jl") end
@time @safetestset "Testing Wish" begin include("TestWish/TestingWish.jl") end
@time @safetestset "Testing cholx and randx" begin include("TestingCholXRandX.jl") end
@time @safetestset "Testing some functions in extfunc" begin include("TestingExtfunc.jl") end
@time @safetestset "Testing ACF function" begin include("testAcf/testAcf.jl") end
@time @safetestset "Testing decVCV function" begin include("testdecVCV.jl") end
@time @safetestset "Testing the fitting of distributions" begin include("testFitDist.jl") end
@time @safetestset "Testing the data manipulation function" begin include("test_data_manipulation.jl") end
@time @safetestset "Testing generate data functions" begin include("test_generate.jl") end
@time @safetestset "Testing PIT estimation" begin include("testPIT.jl") end



