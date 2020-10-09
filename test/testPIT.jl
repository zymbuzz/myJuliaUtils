# testing PIT calculation
using myJuliaUtils
using StatsBase
using Test
using BenchmarkTools


a = randn(10, 10000)
b = zeros(10)
Q = pspit(a, b)
@test length(Q.z) == length(b)
@test sum(Q.zhist.weights) == 10

nbins = rand(1:20)
a = randn(100, 10000)
b = zeros(100)
Q = pspit(a, b, nbins = nbins)
@test length(Q.z) == length(b)
@test length(Q.zhist.weights) == nbins
@test sum(Q.zhist.weights) == 100
# @btime pspit(a,b)

nbins = rand(1:20)
d = rand(100:200)
ndraws = rand(10000:20000)
a = randn(d, ndraws)
b = randn(d)
Q = pspit(a, b, nbins = nbins)
@test length(Q.z) == length(b)
@test length(Q.zhist.weights) == nbins
@test sum(Q.zhist.weights) == d

# Q = pspit(a, b, nbins = nbins, returnZ=true)

# using DelimitedFiles
# AAA=readdlm("/Users/zymantas/Downloads/testPIT.txt", ',');
# pspit(AAA[:,1:end-1], AAA[:,end])