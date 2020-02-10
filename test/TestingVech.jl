using myJuliaUtils
using Test
using LinearAlgebra

randomdim = rand(1:10)
c = Int(randomdim * (randomdim + 1) / 2)
a = randn(c)
b = @inferred vec2sym(a)
@test a == vech(b)
b = vec2ltri(a)
@test a == vech(b)

reps = 1000
q = randn(c, reps)
w = vec2sym(q)
@test w == vec2sym(q', 2)

e = similar(q)
for i = 1:reps
    e[:,i] = vech(w[:,:,i])
end
@test e == q


# testing vec2ltriW1 and vech(,-1)
v = [0.3]
B = [1. 0. ; 0.3 1.]
@test vec2ltriW1(v) == B
@test vec2ltriW1(v, 2) == B
@test vech(B, -1) == v

v = [0.5, 0.3, 0.2]
B = [1. 0. 0.; 0.5 1. 0.; 0.3 0.2 1.]
@test vec2ltriW1(v) == B
@test vec2ltriW1(v, 2) == B
@test vech(B, -1) == v

v = [4., 3., 2., 1., 5., 7.]
B = [1.0  0.0  0.0  0.0;
4.0  1.0  0.0  0.0;
3.0  2.0  1.0  0.0;
1.0  5.0  7.0  1.0]

A = [1.0  0.0  0.0  0.0;
4.0  1.0  0.0  0.0;
3.0  1.0  1.0  0.0;
2.0  5.0  7.0  1.0]
@test vec2ltriW1(v, 2) == B
@test vec2ltriW1(v, 1) != B
@test vec2ltriW1(v, 1) == A
@test vech(A, -1) == v
@test vech(B, -1) != v

randomdim = Int(ceil.(rand() * 20)) + 1
c = Int(randomdim * (randomdim - 1) / 2)
a = randn(c)
b = @inferred vec2ltriW1(a)
@test a == vech(b, -1)

# testing the length of results to be all fine

randomdim = Int(ceil.(rand() * 20)) + 1
c = Int(randomdim * (randomdim - 1) / 2)
a = randn(c)
@test diag(vec2ltriW1(a)) == ones(randomdim)

randomdim = Int(ceil.(rand() * 20))
c = Int(randomdim * (randomdim + 1) / 2)
a = randn(randomdim, randomdim)
@test length(vech(a)) == c

randomdim = Int(ceil.(rand() * 20)) + 1
c = Int(randomdim * (randomdim - 1) / 2)
a = randn(randomdim, randomdim)
@test length(vech(a, -1)) == c


# testing vechRshp
N=rand(1:20)
b=randn(N,N,N);
for i=1:N
    b[:,:,i]=genPSDmat(N)
end

A=@inferred vechRshp(b)
@test A[:,end]==vech(b[:,:,end])
@test A[:,1]==vech(b[:,:,1])
c=rand(1:N)
@test A[:,c]==vech(b[:,:,c])
