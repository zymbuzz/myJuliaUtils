using TVPVARPkg
using Test
using LinearAlgebra
using DelimitedFiles
using Statistics

path = dirname(@__FILE__)
# path = "/Users/zymantas/TVPVARPkg/test/testAcf"

y = readdlm("$path/y.txt", ',');
acf2match = readdlm("$path/acf.txt", ',');
    
@test acf(y[:,1], 0:40) ≈ acf2match[:,1]

@test acf(y, 0:40) ≈ acf2match

for i=1:40
    @test acf(y, i) ≈ acf2match[i+1,:]'
end