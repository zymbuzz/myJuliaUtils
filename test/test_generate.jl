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

        if mod(i,100000)==0
            println("keep")
        end

    end
    lambbdaMean ./= K

    @test isapprox(lambbdaMean, lammbda, atol = 0.01)

end


@testset "testing genMC function" begin
    noStates = round(Int, 2 + rand() * 50)
    S0 = floor(Int, 1 + rand() * noStates)
    PI = eye(noStates)
    T = 100

    @test all(simMC(PI, T, S0=S0) .== S0)

    noStates = round(Int, 2 + rand() * 20)
    S0 = floor(Int, 1 + rand() * noStates)

    A = randn(noStates, noStates)
    A = abs.(A)
    A = A ./ sum(A;dims = 2)
    B = copy(A)
    A[end,:] = [zeros(1, noStates - 1) 1]
    B[1,:] = [1 zeros(1, noStates - 1)]

    F = simMC(A, 1000, S0=S0)
    @test F[end] == noStates
    @test F[1] == S0

    F = simMC(B, 1000, S0=S0)
    @test F[end] == 1
    @test F[1] == S0
end