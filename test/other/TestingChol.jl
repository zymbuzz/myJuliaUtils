
using myJuliaUtils
using Test
using LinearAlgebra
using Statistics

m = 0
reps = 5
for i = 1:reps
    K = round(Int, 1 + rand() * 9);
    uV = randn(K, K);
    uV = uV'uV
    epsilon = cholesky(uV).L * randn(K, 50000000);
    global m += cov(epsilon') ≈ uV
end
@test m == reps