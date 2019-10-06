using Test
using myJuliaUtils

T = 100;
K = rand(1:20);
x = abs.(randn(T, K))
@testset "testing dataTransf" begin

    a = dataTransf(ones(T), 3)
    @test a.y == zeros(T)
    @test a.n_loss == 0

    @test dataTransf(x, 1) == dataTransf(x, ones(Int, K))


    transfId = rand(1:6)
    (a1, b1) = dataTransf(x, transfId)
    (a2, b2) = dataTransf(x, transfId * ones(Int, K))

    @test b1 == b2
    @test a1[b1 + 1:end,:] == a2[b2 + 1:end,:]
end