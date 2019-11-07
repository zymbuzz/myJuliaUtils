function lag0(x::AbstractVector, p::Int = 1)
    T = length(x);
    return [zeros(eltype(x), p); x[1:T - p]]
end

function lag0(x::AbstractArray, p::Int = 1)
    T, n = size(x);
    return [zeros(eltype(x), p, n); x[1:T - p,:]]
end

function eye(m::Int, n::Int = m)
    return Matrix{Float64}(I, m, n)
end

maxabs(x) = maximum(abs, x)

function vech(A::AbstractMatrix, k::Int = 0)
    # return A[tril!(trues(size(A)), 0)]
    return A[tril(trues(size(A)), k)]
end

function vech(A::AbstractArray{T,3}, k::Int = 0) where {T}
    B = vech(A[:,:,1], k)
    a = size(A, 3)

    if a > 1
        for i = 2:a
            append!(B, vech(A[:,:,i], k))
        end
    end

    return B
end

# this one is type unstable
function vec2ltri2(v::AbstractVector)
    n = length(v)
    s = round(Int, (sqrt(8n + 1) - 1) / 2)
    @assert s * (s + 1) / 2 == n "length of vector is not triangular"
    k = 0; [ i >= j ? (k += 1; v[k]) : zero(eltype(v)) for i = 1:s, j = 1:s ]
end

function vec2ltri(v::AbstractVector)
    n = length(v)
    s = round(Int, (sqrt(8n + 1) - 1) / 2)
    @assert s * (s + 1) / 2 == n "length of vector is not triangular"
    A = zeros(eltype(v), s, s)
    k = 0
    for j = 1:s, i = j:s
        k += 1
        A[i,j] = v[k]
    end
    return A
end

function vec2ltriW1(v::AbstractVector{T}, dm::Int = 1)::Array{T,2} where {T}
    n = length(v)
    s = round(Int, (sqrt(8n + 1) + 1) / 2)
    @assert s * (s - 1) / 2 == n "length of vector is not triangular below diagonal "
    A = Matrix{T}(I, s, s)
    k = 0
    if dm == 1
        for j = 1:s, i = j:s
            if i > j
                k += 1  
                A[i,j] = v[k]
            end
        end
    elseif dm == 2
        for i = 1:s, j = 1:i
            if i > j
                k += 1  
                A[i,j] = v[k]
            end
        end
    end
    return A
end

function vec2sym(v::AbstractVector)
    return Symmetric(vec2ltri(v), :L)
end

function vec2sym(v::AbstractMatrix, dm::Int = 1)
    if dm == 2
        v = v'
    end

    a, b = size(v)
    s = round(Int, (sqrt(8a + 1) - 1) / 2)

    A = similar(v, s, s, b)
    for i = 1:b
        A[:,:,i] = vec2sym(v[:,i])
    end
    return A
end

function companionf(beta::AbstractMatrix, n::Int, p::Int, c::Bool = false)
    return [beta[:, 1:end - c];eye(n * p - n) zeros(n * p - n, n)]
end

function companionf(beta::AbstractVector, n::Int, p::Int, c::Bool = false)
    @assert n == 1 throw(ArgumentError("provided beta is univariate"))
    if c || p > 1
        return [reshape(beta[1:end - c], 1, :); eye(p - 1) zeros(p - 1)]
    else
        return Matrix{Float64}(reshape(beta, 1, 1))
    end  
end

function stabcheck(beta::AbstractVecOrMat, n::Int, p::Int, c::Bool = false)::Bool
    # beta is (Nx(N*L+c))
    FF = companionf(beta, n, p, c)
    return maximum(abs, eigvals(FF)) >= 1.
end

function stabcheckC(beta::AbstractVecOrMat, n::Int, p::Int, c::Bool = false)::Bool
    # beta is (Nx(N*L+c))
    # stability check which ignores complex numbers
    FF = companionf(beta, n, p, c)
    D = eigvals(FF)
    if isreal(D)
        S = maximum(abs, D) >= 1.
    else
        S = true
    end
    return S
end

function stabcheck(beta::Number, n::Int = 1, p::Int = 1, c::Bool = false)::Bool
    return abs(beta) >= one(eltype(beta))
end
stabcheckC(beta::Number, n::Int = 1, p::Int = 1, c::Bool = false) = stabcheck(beta, n, p, c)

function preparexy(data::AbstractVector, p::Int, c::Bool = false)
    T = length(data);
    Y = data[p + 1:end];
    X = similar(data, T, p + c);
    for i = 1:p
        X[:,(i - 1) + 1:i] = lag0(data, i);
    end
    if c
        X[:,end] = ones(eltype(data), T);
    end
    X = X[p + 1:end,:]
    return Y, X
end

function preparexy(data::AbstractMatrix, p::Int, c::Bool = false)
    T, n = size(data);
    Y = data[p + 1:end,:];
    X = similar(data, T, n * p + c);
    for i = 1:p
        X[:,n * (i - 1) + 1:n * i] = lag0(data, i);
    end
    if c
        X[:,end] = ones(eltype(data), T);
    end
    X = X[p + 1:end,:]
    return Y, X
end

function sumsqr(a::AbstractArray)
    return mapreduce(x->x^2,+,a)
end

# function wish(h, n::Int)
#     A = cholesky!(h).L * randn(size(h, 1), n)
#     A = A * A'
#     return Symmetric(A)
# end

function wish(h::AbstractArray{S,2}, n::Int)   where {S <: AbstractFloat}
    A = randn(n, size(h, 1)) * cholesky(h).U
    return A'A
end

function wish(h::AbstractArray{S,1}, n::Int)   where {S <: AbstractFloat}
    A = randn(n) .* sqrt.(h)
    return dot(A', A)
end

function iwish(h::AbstractMatrix, n::Int)
    B = randn(n, size(h, 1))
    A = cholesky(h).L / qr!(B).R
    return A * A'
end

function iwish(h::AbstractVector, n::Int)
    B = randn(n, 1)
    A = sqrt.(h) / qr!(B).R
    return dot(A, A')
end

function genPSDmat(K::Int = rand(2:9))
    A = randn(K, K)
    return A'A
end

function genPSDmatStrict(K::Int = rand(2:9))
    A = randn(K, K)
    A = A'A
    vals, vecs = eigen(A)
    vals[1] = 0.
    return vecs * Diagonal(vals) * vecs'
end

function genPDmat(K::Int = rand(2:9))
    A = zeros(K, K)
    while isposdef(A) == false
        A = randn(K, K)
        A = A'A
    end
    return A
end

function cholPSD(A::AbstractMatrix{T})::Array{T,2}  where {T <: AbstractFloat}
    B = copy(Symmetric(A));
    if isposdef(B)

        # B = cholesky(B)
        # return convert(Array{Float64,2}, B.U)

        return cholesky!(B).U

    else
        # A = sqrt(A)
        # # D,U = eig(A);

        # A = eigen(A);
        
        # tol = length(A.values) * eps(maximum(A.values));
        # # D1 = A.values[A.values .> tol];
        # A = diagm(sqrt.(A.vectors)) * A.vectors[:, A.values .> tol]';

        # M = eigen(B)
        # D = M.values
        # U = M.vectors

        D, U = eigen(B)
        tol = length(D) * eps(maximum(D));
        # D1 = D[D .> tol];
        # # A = diagm(0 => sqrt.(D1)) * U[:, D .> tol]';
        # C = Diagonal(sqrt.(D1)) * U[:, D .> tol]';

        # return Diagonal(sqrt.(filter(z -> z > tol, D))) * U[:, D .> tol]'
        cord = findall(z->z > tol, D)
        # return Diagonal(sqrt.(filter(z -> z > tol, D))) * U[:, findall(z -> z > tol, D)]'
        return Diagonal(sqrt.(D[cord])) * U[:, cord]'

    end
end

function randnPSD(mu::AbstractMatrix{T}, sigma::AbstractMatrix{T})::Array{T}  where {T <: AbstractFloat}
    # input:
    # mu is row vector
    # 
    # output:
    # draw is a row vector
    sigmaN = cholPSD(sigma)
    draw = mu + randn(1, size(sigmaN, 1)) * sigmaN 
    return draw
end

function randnPSD(mu::AbstractVector{T}, sigma::AbstractMatrix{T})::Array{T}  where {T <: AbstractFloat}
    # input:
    # mu is row vector
    # 
    # output:
    # draw is a row vector
    sigmaN = cholPSD(sigma)
    draw = mu' + randn(1, size(sigmaN, 1)) * sigmaN 
    return draw
end

function mcSS(PI::Matrix)
    a = eigen(Matrix(PI'))
    b = a.vectors[:,a.values .≈ 1.]
    return vec(b) ./ sum(b)
end

function mcSS2(PI::Matrix, tol::Float64 = eps())
    l= size(PI, 1)
    q = ones(l)./l
    q1 = similar(q)
    a = 2.
    while a > tol
        q1 .= PI'*q
        a = maximum(abs, q1 .- q)
        q .= q1
    end
    return q
end

function quantileArr(v::AbstractArray{T,2}, p, dms::Int = 1) where {T}
    a, b = size(v)
    w = length(p)
    p = vec(p)
    if dms == 1
        A = similar(v, w, b)
        for i = 1:b
            A[:,i] = quantile(v[:,i], p)
        end
    elseif dms == 2
        A = similar(v, a, w)
        for i = 1:a
            A[i,:] = quantile(v[i,:], p)
        end
    end
    return A
end

function quantileArr(v::AbstractArray{T,3}, p, dms::Int = 1) where {T}
    a, b, c = size(v)
    w = length(p)
    p = vec(p)

    if dms == 1
        A = similar(v, w, b, c)
        for i = 1:b, j = 1:c
            A[:,i,j] = quantile(v[:,i,j], p)
        end
    elseif dms == 2
        A = similar(v, a, w, c)
        for i = 1:a, j = 1:c
            A[i,:,j] = quantile(v[i,:,j], p)
        end
    elseif dms == 3
        A = similar(v, a, b, w)
        for i = 1:a, j = 1:b
            A[i,j,:] = quantile(v[i,j,:], p)
        end
    end

    return A
end

ismyapprox(a, b, dgts::Int) = isapprox(a, b, atol = dgts)

function normpdf(x::Float64, mu::Float64 = 0.0, sigma::Float64 = 1.0)
    return exp(-0.5 * ((x - mu) ./ sigma).^2) ./ (sqrt(2. * pi) .* sigma)
end
      
function acf(y::AbstractVector, q::Int)
# function assumes that there are no missings

    z = y .- mean(y)
    T = length(y)

    ac0 = dot(z, z) ./ T
    ac = dot(z[1:end - q], z[q + 1:end]) ./ T
    ac /= ac0; # Normalize to correllation

    return ac
end

function acf(y::AbstractVector, q::AbstractRange = 1:10)
    # function assumes that there are no missings
    
    z = y .- mean(y)
    T = length(y)

    ac0 = dot(z, z) ./ T
    ac = Array{Float64}(undef, 0);
    for j = q
        push!(ac, dot(z[1:end - j], z[j + 1:end]) ./ T)
    end
    ac ./= ac0; # Normalize to correllation

    return ac
end

function acf(y::AbstractMatrix, q = 1:10)

    T, n = size(y)

    ac = Array{Float64}(undef, 0);

    for i = 1:n
        append!(ac, acf(y[:,i], q))
    end

    ac = reshape(ac, length(q), n)

    return ac
end

function decVCV(VCV::AbstractVector, getAinv::Bool = false)
    Omega = vec2sym(VCV)
    Omega = cholesky(Omega).L
    sigma = diag(Omega)
    A = Omega / Diagonal(sigma)
    if getAinv
        A = inv(A)
    end
    A = vech(A, -1)

    return A, sigma
end

function decVCV(VCV::Matrix{Float64}, getAinv::Bool = false)

    a, T = size(VCV)

    s = round(Int, (sqrt(8a + 1) - 1) / 2)
    s * (s + 1) / 2 == a || error("fromVCV: length of vector is not triangular")

    sigma = Array{Float64}(undef, s, T)
    A = Array{Float64}(undef, round(Int, s * (s - 1) / 2), T)

    for i = 1:T
        A[:,i], sigma[:,i] = decVCV(VCV[:,i], getAinv)
    end

    return A, sigma
end

function ols1(y::AbstractVector, x::AbstractMatrix)
    nobs, nvar = size(x)
    bhat = x \ y
    yhat = x * bhat
    res = y - yhat
    sig2hat = (res'res) / (nobs - nvar)
    sigbhat = sig2hat .* eye(nvar) / (x' * x)
    R2 = var(yhat) / var(y)
    
    return (nobs = nobs, nvar = nvar, bhat = bhat, yhat = yhat, res = res, sig2hat = sig2hat, sigbhat = sigbhat, R2 = R2)
end

function VARols1(data::AbstractVecOrMat{Float64}, cnst::Bool, P::Int)
    y, x = preparexy(data, P, cnst)
    nobs, nvar = size(x)
    bhat = x \ y
    yhat = x * bhat
    res = y - yhat
    sig2hat = (res'res) / (nobs - nvar)
    
    return (nobs = nobs, nvar = nvar, bhat = bhat, yhat = yhat, res = res, sig2hat = sig2hat)
end

function nanmean(x::AbstractArray{T}) where T <: AbstractFloat
    sum = zero(eltype(x))
    count = 0
    for i in x
        if !isnan(i)
            sum += i
            count += 1
        end
    end
    return sum / count
end

nanmean(x::AbstractArray, y::Int) = mapslices(nanmean, x; dims = y)

function inbetween(a::Matrix,b::AbstractArray{T,3})  where {T}
    return (view(b,:,:,1) .<= a) .& (a.<= view(b,:,:,2))
end
