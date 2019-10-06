using myJuliaUtils
using Test
using LinearAlgebra
using DelimitedFiles

# Testing both lag functions
@testset "testing lag functions" begin
    # if vector
    T = 100
    P = rand(1:T - 1)
    a = collect(1.0:1.0:T)
    b = lag0(a, P)
    @test b[1:P] == zeros(P)
    @test T - P == b[end]

    # if matrix
    a = reshape(a, T, 1)
    P = round(Int, 1 + rand() * T - 1)
    b = lag0(a, P)
    @test b[1:P,1] == zeros(P)
    @test T - P == b[end,1]

    a = [ones(T, 1) a]
    P = round(Int, 1 + rand() * T - 1)
    b = lag0(a, P)
    @test b[1:P,2] == zeros(P)
    @test T - P == b[end,2]
end

@testset "testing function preparexy()" begin
    T = 100
    P = rand(1:T - 1)
    a = randn(T)
    cnst = rand() > 0.5
    Y, X = preparexy(a, P, cnst)
    # for univariate
    i = rand(1:T - P)
    z = reverse(a[i:i + P - 1])
    if cnst
        push!(z, 1.)
    end
    @test X[i,:] == z
    @test Y[i] == a[i + P]

    # for multivariate
    N = rand(1:20)
    a = randn(T, N)
    Y, X = preparexy(a, P, cnst)
    nT, nX = size(X)

    i = rand(1:T - P)
    z = Vector(vec(reverse(a[i:i + P - 1,:], dims = 1)'))
    if cnst
        push!(z, 1.)
    end
    @test nT == T - P
    @test nX == P * N + cnst
    @test X[i,:] == z
    @test Y[i,:] == a[i + P,:]
end


# testing companion form function
@testset "testing companion form functions" begin
    betta = [2 3;4 5]
    n = 2
    p = 1
    c = false
    @test betta == companionf(betta, n, p, c)

    bettaC = [betta ones(n, 1)]
    c = true
    @test betta == companionf(bettaC, n, p, c)

    c = false
    p = 2
    betta = [2 3 4 5;6 7 8 9]
    betta_comp = [betta;1 0 0 0; 0 1 0 0]
    @test betta_comp == companionf(betta, n, p, c)

    bettaC = [betta ones(n, 1)]
    c = true
    @test betta_comp == companionf(bettaC, n, p, c)

    betta = [1.]
    n = 1
    p = 1
    c = false
    @test reshape(betta, 1, 1) == companionf(betta, n, p, c)

    n = 2
    p = 1
    c = false
    @test_throws ArgumentError companionf(betta, n, p, c)

    n = 1
    p = 2
    c = false
    @test_throws ArgumentError companionf(betta, n, p, c)

    n = 1
    p = 1
    c = true
    @test_throws ArgumentError companionf(betta, n, p, c)
end

# Testing stability function
@testset "testing stability functions" begin
    N = round(Int, 1 + rand() * 100)
# A = randn(N, N)
# A = A'A
    A = genPDmat(N)
    S, Q = eigen(A)
    S[end] = 1.0
    A = Q * Diagonal(S) * Q'
    @test stabcheck(A, N, 1, false)
    @test stabcheckC(A, N, 1, false)

    a = findall(x->x >= 1.0, S)
    b = length(a)
    S[a] = ones(b) - rand(b)
    A = Q * Diagonal(S) * Q'
    @test stabcheck(A, N, 1, false) == false
    @test stabcheckC(A, N, 1, false) == false

    A = [   
    0.516090886728283	0.0724233649560897	0.0816329971852510	0.143587089009963	-0.00458107252041236
    0.0724233649560897	0.844836851975328	0.176990845416178	-0.0777482352998716	-0.0921908525982972
    0.0816329971852510	0.176990845416178	0.671194897897464	-0.00975809606339381	-0.0331333053292128
    0.143587089009963	-0.0777482352998717	-0.00975809606339381	0.797277510951764	-0.218043826201017
    -0.00458107252041236	-0.0921908525982972	-0.0331333053292128	-0.218043826201017	0.554361340591299
];

    @test stabcheck(A, size(A, 1), 1, false)
    @test stabcheckC(A, size(A, 1), 1, false)

    A = [
    0.500514390972542	-0.0693316507737557	0.103250184760922	0.157734480886660	0.0498798880488583
    -0.0693316507737557	0.586863186538567	0.00546111896077394	-0.0714652123701029	0.153497763972517
    0.103250184760922	0.00546111896077393	0.487346580847464	-0.0744321935787715	-0.0521727141915581
    0.157734480886660	-0.0714652123701029	-0.0744321935787715	0.566439153014101	0.0799567015086364
    0.0498798880488583	0.153497763972517	-0.0521727141915581	0.0799567015086364	0.603872518537494
];
    @test stabcheck(A, size(A, 1), 1, false) == false
    @test stabcheckC(A, size(A, 1), 1, false) == false


# AR(2) stability check
    Nvar = 1;
    LAGS = 2;
    c = false;
    reps = 50000
    BETTAt = ones(LAGS, reps);
    for i = 1:reps
        chk = 1
        while chk == true
            global BETTAt[:,i] = randn(LAGS, 1)
            chk = stabcheck(BETTAt[:,i]', Nvar, LAGS, c)
        end
    end

    @test all([BETTAt[2,:] .< BETTAt[1,:] .+ 1. BETTAt[2,:] .< -BETTAt[1,:] .+ 1. BETTAt[2,:] .> -1. ])

# some other small tests

    @test stabcheck(1., 1, 1, false) == true
    @test stabcheck([1.], 1, 1, false) == true
    @test stabcheck(ones(1, 1), 1, 1, false) == true
    @test stabcheck([1. randn()], 1, 1, true) == true
    @test stabcheckC(1., 1, 1, false) == true
    @test stabcheckC([1.], 1, 1, false) == true
    @test stabcheckC([1. randn()], 1, 1, true) == true
    @test stabcheck(1, 1, 1, false) == true
    @test stabcheck([1 randn()], 1, 1, true) == true
    @test stabcheckC(1, 1, 1, false) == true
    @test stabcheckC([1 randn()], 1, 1, true) == true
end

# testing simMC [to add more tests]

noStates = round(Int, 2 + rand() * 50)
S0 = floor(Int, 1 + rand() * noStates)
PI = eye(noStates)
T = 100

@test all(simMC(PI, T, S0) .== S0)

noStates = round(Int, 2 + rand() * 20)
S0 = floor(Int, 1 + rand() * noStates)

A = randn(noStates, noStates)
A = abs.(A)
A = A ./ sum(A;dims = 2)
B = copy(A)
A[end,:] = [zeros(1, noStates - 1) 1]
B[1,:] = [1 zeros(1, noStates - 1)]

F = simMC(A, 1000, S0)
@test F[end] == noStates
@test F[1] == S0

F = simMC(B, 1000, S0)
@test F[end] == 1
@test F[1] == S0


# test for normal pdf
x = [-2.,-1.,0.,1.,2.]
y1 = [0.0539909665131881, 0.241970724519143, 0.398942280401433, 0.241970724519143, 0.0539909665131881]

mu = 2.
sigma = 1.
y2 = [0.000133830225764885, 0.00443184841193801, 0.0539909665131881, 0.241970724519143, 0.398942280401433]

@testset "testing normal pdf function" begin
    @test normpdf.(x) ≈ y1
    @test normpdf.(x, mu, sigma) ≈ y2
end



# testing ols function
path = dirname(@__FILE__)
# path = "/Users/zymantas/TVPVARPkg/test/"

data = readdlm("$path/TestOls/data2test.txt", ',');
B2check = readdlm("$path/TestOls/B2check.txt", ',');
stats2check = readdlm("$path/TestOls/stats2check.txt", ',');
covb2check = readdlm("$path/TestOls/covb2check.txt", ',');

y = data[:,1]
x = [ones(length(y)) data[:,2:end]]
r = ols1(y, x)

@testset "testing OLS function" begin
    @test r.bhat ≈ B2check
    @test r.sigbhat ≈ covb2check
    @test r.R2[1] ≈ stats2check[1]
    @test r.sig2hat[1] ≈ stats2check[2]
end

# testing quantileArr function
a = randn(rand(1:20), rand(1:20), rand(1:20))

@testset "testing quantileArr function" begin
    @test quantileArr(a, [0.5], 1) ≈ median(a, dims = 1)
    @test quantileArr(a, [0.5], 2) ≈ median(a, dims = 2)
    @test quantileArr(a, [0.5], 3) ≈ median(a, dims = 3)
end

@testset "testing nan mean" begin
    @test nanmean([1., 2., NaN]) == 1.5
    @test nanmean([1. 2.; NaN 3.]) == 2.0
    @test nanmean([1. 2.; NaN 3.], 2) ≈ [1.5; 3. ]
    @test nanmean([1. 2.; NaN 3.], 1) ≈ [1.0 2.5 ]
    a = [1. 2.; NaN 3.]
    a = cat(a, [1. NaN; 7. 3.], dims = 3)
    @test nanmean(a, 3) ≈ [1. 2.; 7. 3.]


    a = randn(4, 3, 5)
    @test mean(a;dims = 1) ≈ nanmean(a, 1)
    @test mean(a;dims = 2) ≈ nanmean(a, 2)
    @test mean(a;dims = 3) ≈ nanmean(a, 3)

end


@testset "testing inbetween" begin
    b = randn(2, 2, 2)
    b[:,:,1] = [1 2; 3 4]
    b[:,:,2] = b[:,:,1] .+ 2
    a = [3 5; 9 5]
    @test inbetween(a, b) == [true false; false true]

    T = rand(1:1000)
    n = rand(1:1000)

    a = randn(T, n)
    b = randn(T, n, 2)
    int = abs(randn())
    b[:,:,2] = b[:,:,1] .+ int
    a = b[:,:,1] .+ int ./ 2

    @test all(inbetween(a, b)) 

end


