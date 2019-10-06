using TVPVARPkg
using Test
using LinearAlgebra

N = rand(1:100)
A = genPDmat(N)
cholA = cholPSD(A)
@test cholA == cholesky(A).U
@test cholA'cholA ≈ A

# test1=0
# test2=0
# NoTimes=20
# for i=1:NoTimes
# N = round(Int, 1 + rand() * 100)
# A = randn(N, N)
# A = A'A
# cholA = cholPSD(A)
# test1+=cholA == cholesky(A).U
# test2+=cholA'cholA ≈ A
# end
# @test test1==NoTimes
# @test test2==NoTimes

N = rand(2:100)
K = rand(1:N - 1)
A = genPDmat(N)
S, Q = eigen(A)
# S[1:K] = ones(K) * eps(0.0)
S[1:K] = zeros(K)
# A = Q * Diagonal(S) * inv(Q)
A = Q * Diagonal(S)
A[:,1:K] = zeros(N, K)
A .= A / Q
cholA = cholPSD(A)
if K > 1
    @test size(cholA) == (N - K, N)
else 
    @test all(size(cholA) .>= (N - K, N)) # this is done to avoid approximation error
end

@test cholA'cholA ≈ A


# Nseed=round(Int, 1 + rand() * 5000)
# srand(Nseed)
# a=randn(1, N) * cholx(A)
# srand(Nseed)
# a=randn(1, N) * cholx(A)