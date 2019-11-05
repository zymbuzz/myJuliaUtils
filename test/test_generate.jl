using myJuliaUtils
using Test
using LinearAlgebra
using DelimitedFiles
using Statistics

@testset "testing genLPM function - Monte Carlo test" begin
    T = 10000
    K = 200000
    X = [ones(T) randn(T)];
    lammbda = [rand(),randn() / 10.]

    lambbdaMean = 0. .* lammbda
    y = BitArray(undef, T)
    @time for i = 1:K
        y .= genLPM(lammbda, T, X)
        lambbdaMean .+= X \ y

        if mod(K,10000)==0
            println("keep")
        end

    end
    lambbdaMean ./= K

    @test isapprox(lambbdaMean, lammbda, atol = 0.01)

end