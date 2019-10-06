using TVPVARPkg
using Test
using LinearAlgebra
using Random

# testing some helper function

# one matrix case (qq is a vector)
n = round(Int, 2 + rand() * 100)
qq = genPDmat(n)
qq = vech(qq)

A, sigma = TVPVARPkg.decVCV(qq)

@test vec2sym(qq) ≈ vec2ltriW1(A) * Diagonal(sigma).^2 * vec2ltriW1(A)'

# many matrix case (qq is a matrix)

n = 5
T = 100
c = round(Int, n * (n + 1) / 2)
qq = Array{Float64}(undef, c, T)

for i = 1:T
    q = genPDmat(n)
    qq[:,i] = vech(q)
end

A, sigma = TVPVARPkg.decVCV(qq)

for i = 1:T
    @test vec2sym(qq[:,i]) ≈ vec2ltriW1(A[:,i]) * Diagonal(sigma[:,i]).^2 * vec2ltriW1(A[:,i])'
end

A, sigma = TVPVARPkg.decVCV(qq,true)

for i = 1:T
    @test vec2sym(qq[:,i]) != vec2ltriW1(A[:,i]) * Diagonal(sigma[:,i]).^2 * vec2ltriW1(A[:,i])'
end
