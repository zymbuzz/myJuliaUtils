using myJuliaUtils
using Test
using LinearAlgebra
using Random
using Distributions
# using BenchmarkTools
using DelimitedFiles

# comparing the procedure to matlab

path = dirname(@__FILE__)
# path = "/Users/zymantas/TVPVARPkg/test/TestWish"


T = readdlm("$path/T.txt", ',');
T = Int(T[1])
if T > 80
    error("matlab uses a different procedure for high degrees of freedom")
end
draw = readdlm("$path/draw.txt", ',');
aa = readdlm("$path/WW.txt", ',');
MatlabInvWish = readdlm("$path/matlabinvwish.txt", ',');
MatlabWish = readdlm("$path/matlabwish.txt", ',');

A = draw * cholesky(aa).U
JuliaWish = A'A


B = copy(draw)
A = cholesky(aa).L / qr!(B).R
JuliaInvWish = A * A'

# C2 = inv(Symmetric(aa));
# C3= draw * cholesky(C2).U
# C3=C3'C3
# D2 = inv(Symmetric(C3))
# @test D2 ≈ MatlabInvWish

@testset "Compared to Matlab" begin
    @test JuliaWish ≈ MatlabWish
    @test JuliaInvWish ≈ MatlabInvWish
end


# testing wishart and inverse wishart for univariate compared to multivariate
WW = randn(1, 1)
WW = WW'WW
T = rand(1:1000)

SomeSeed = rand(1:10000)

Random.seed!(SomeSeed);
C = wish(WW, T) 

Random.seed!(SomeSeed);
C1 = wish(vec(WW), T)



Random.seed!(SomeSeed);
D = iwish(WW, T) 

Random.seed!(SomeSeed);
D1 = iwish(vec(WW), T)



@testset "matching univariate to multivariate" begin
    @test C ≈ [C1]
    @test D ≈ [D1]
end

# testing wishart and inverse wishart for univariate
WW = randn(1, 1)
WW = WW'WW
T = rand(1:1000)

Random.seed!(SomeSeed);
A = rand(Wishart(T, WW))
Random.seed!(SomeSeed);
C = wish(WW, T)

Random.seed!(SomeSeed);
A = rand(InverseWishart(T, WW))

Random.seed!(SomeSeed);
B = iwish(WW, T)

Random.seed!(SomeSeed);
C2 = inv(WW);
C3 = wish(C2, T);
C = inv(C3)


# for univariate case where WW is a vector
WW = [1 + rand() * 10000]
T = rand(1:1000)

SomeSeed = rand(1:10000)

Random.seed!(SomeSeed);
A_v = rand(InverseWishart(T, reshape(WW, 1, 1)))

Random.seed!(SomeSeed);
B_v = iwish(WW, T)

Random.seed!(SomeSeed);
C2 = inv.(WW);
C3 = wish(C2, T);
C_v = inv.(C3)

@testset "univariate" begin
    # @test A ≈ B
    @test B ≈ C

    # @test A_v ≈ B_v
    @test B_v ≈ C_v
end



# # testing for multivariate wishart

# K = round(Int, 1 + rand() * 6)

# WW = randn(K, K)
# WW = WW'WW
# T = round(Int, 1 + rand() * 1000)

# # testing the equality between wish and iwish
# Random.seed!(SomeSeed);
# A = rand(InverseWishart(T, WW))
# Random.seed!(SomeSeed);
# B = iwish(WW, T)
# Random.seed!(SomeSeed);
# C2 = inv(Symmetric(WW));
# C3 = wish(C2, T);
# C = inv(C3)

# @testset "multivariate" begin
#     @test_broken A ≈ B
#     @test_broken B ≈ C
# end



# using Law of large numbers

reps = 2000000
K = rand(1:5)

WW = randn(K, K)
WW = WW'WW
T = rand(K+1:300)

SomeSeed = rand(1:10000)
Random.seed!(SomeSeed)
a1 = wish(WW, T)

Random.seed!(SomeSeed);
a2 = rand(Wishart(T, WW))

cont1 = zeros(K, K);
cont2 = copy(cont1);
cont3 = copy(cont1);
cont4 = copy(cont1);

wD=Wishart(T, WW)
iwD=InverseWishart(T, WW)

for i = 1:reps
    global cont1 .+= wish(WW, T);
    global cont2 .+= rand(wD);

    global cont3 .+= iwish(WW, T);
    global cont4 .+= rand(iwD);
end
cont1 ./= reps;
cont2 ./= reps;
cont3 ./= reps;
cont4 ./= reps;

@testset "based on LLN" begin

    # @test maximum(abs.(round.(cont1 - cont2; digits=0)))==0.0
    @test maximum(abs.(round.(cont1 - mean(Wishart(T, WW)); digits=0)))==0.0
    @test maximum(abs.(round.(mean(Wishart(T, WW)) - cont2; digits=0)))==0.0

    # @test maximum(abs.(round.(cont3 - cont4; digits=4)))==0.0
    @test maximum(abs.(round.(cont3 - mean(InverseWishart(T, WW)); digits=3)))==0.0
    @test maximum(abs.(round.(mean(InverseWishart(T, WW)) - cont4; digits=3)))==0.0

end