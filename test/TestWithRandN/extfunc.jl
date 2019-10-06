function lag0(x::Array{Float64}, p::Int64=1)
    T = length(x);
    return [zeros(p, 1); x[1:T - p]]
end

function lag0(x::Array{Float64,2}, p::Int64=1)
    T, n = size(x);
    return [zeros(p, n);x[1:T - p,:]]
end

function eye(m::Int64, n::Int64=m)
    return Matrix{Float64}(I, m, n)
end

function choleskyss(A)
    return cholesky(Symmetric(A))
end

function invss(A)
    return inv(Symmetric(A))
end

function stabcheck(beta, n::Int64, p::Int64, c::Int64=1)
    # beta is (Nx(N*L+c))
    FF = [beta[:, 1:end - c];eye(n * p - n) zeros(n * p - n, n)]
    return eigmax(FF) > 1
end

function stabcheckC(beta, n::Int64, p::Int64, c::Int64=1)
    # beta is (Nx(N*L+c))
    # stability check which ignores complex numbers
    FF = [beta[:, 1:end - c];eye(n * p - n) zeros(n * p - n, n)]
    D = eigvals(FF)
    if isreal(D)
        S = maximum(D) > 1
    else
        S = false
    end
    return S
end

function preparexy(data::Array{Float64,1}, p::Int64, c::Int64=1)
    T = length(data);
    Y = data[p + 1:end];
    X = zeros(T, p + c);
    for i = 1:p
        X[:,(i - 1) + 1:i] = lag0(data, i);
    end
    if c == 1
        X[:,end] = ones(T);
    end
    X = X[p + 1:end,:]
    return Y, X
end

function preparexy(data::Array{Float64,2}, p::Int64, c::Int64=1)
    T, n = size(data);
    Y = data[p + 1:end,:];
    X = zeros(T, n * p + c);
    for i = 1:p
        X[:,n * (i - 1) + 1:n * i] = lag0(data, i);
    end
    if c == 1
        X[:,end] = ones(T);
    end
    X = X[p + 1:end,:]
    return Y, X
end

function wish(h, n::Int64)
    h = cholesky(h)
    A = h.L * randnNEW(size(h, 1), n)
    A = A * A'
    A = Symmetric(A)
    return A
end

function cholx(A)
    B = copy(Symmetric(A));
    if isposdef(B)

        B = cholesky((B + B') / 2)
        return convert(Array{Float64,2}, B.U)

        # this one is more efficient however there is some issue with multiplication 
        # return cholesky!(B)

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
        cord = findall(z -> z > tol, D)
        # return Diagonal(sqrt.(filter(z -> z > tol, D))) * U[:, findall(z -> z > tol, D)]'
        return Diagonal(sqrt.(D[cord])) * U[:, cord]'

    end
end

function randnx(mu::Union{Array{Float64,1},Array{Float64,2},Adjoint{Float64,Array{Float64,1}},Adjoint{Float64,Array{Float64,2}}}, sigma)::Array{Float64}
    # input:
    # mu is row vector
    # 
    # output:
    # draw is a row vector
    sigmaN = cholx(sigma)
    draw = mu + randnNEW(1, size(sigmaN, 1)) * sigmaN 
    return draw
end

function carterkohnMLTV(y::Array{Float64,2}, Z::Array{Float64,2}, Ht::Union{Symmetric{Float64,Array{Float64,2}},Array{Float64,2}}, Qt::Union{Symmetric{Float64,Array{Float64,2}},Array{Float64,2}}, K::Int64, N::Int64, T::Int64, b0, V0::Union{Symmetric{Float64,Array{Float64,2}},Array{Float64,2}}) 

    # y_t=x_t*b_t+e_t where e_t follows (0,H)
    # b_{t}=b_{t-1}+v_t where v_t follows (0,Q)

    # INPUTS:
    # y     - endogeneous variable (NxT)
    # Z     - exogeneous variable (T*N x NP*N) adjusted s.t. Z = kron(x,eye(n));
    # Ht    - VCV of measurement error
    # Qt
    # b0    - prior mean (Kx1)
    # b0V   - prior variance (KxK)
    # K     - number of parameters
    # T     - number of time periods
    # N     - number of variables

    # OUTPUT:
    # bdraw - draws of unobservables (KxT)

    # Variables:
    # Vp    - predicted variance of unobservable (KxK)
    # bp    - predicted value of unobservable (Kx1)
    # Vtt   - updated variance of unobservable (KxK)
    # btt   - updated value of unobservable (Kx1)

    # Notes
    # here I am assuming the random walk assumption, s.T. the matrix F (nelson) is identity, this suggest that predicted value at t is equivalent to updated at t-1 (or btt=bp)

    bp = copy(b0);
    Vp = copy(V0);
    bt = zeros(K, T);
    Vt = zeros(K^2, T);
    #log_lik = 0;
    R = copy(Ht);
    for i = 1:T
        # R = Ht((i-1)*p+1:i*p,:); 
    
        X = Z[(i - 1) * N + 1:i * N,:]; # (N x NP*N)
        cfe = y[:,i] - X * bp;   #conditional forecast error
        f = X * Vp * X' + R;    # variance of the conditional forecast error (NxN)
    
        invf = invss(f);
        # invf = inv(f);
    
        #log_lik = log_lik + log(det(f)) + cfe'*invf*cfe;
    
        #updating
        btt = bp + Vp * X' * invf * cfe;
        Vtt = Vp - Vp * X' * invf * X * Vp;
    
        # predicting for the next iteration
        Vp = Vtt + Qt;
        bp = copy(btt);
    
        bt[:,i] = copy(btt);
        Vt[:,i] = reshape(Vtt, K^2, 1);
    end

    Vtt = reshape(Vt[:,T], K, K);

    # bdraw = zeros(K, T);
    bdraw = Matrix{Float64}(undef, K, T)
    # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
    bdraw[:,T] = randnx(bt[:,T]', Symmetric(Vtt))';
    
    # Backward recurssions
    for i = T - 1:-1:1
        # println(i)
        bf = copy(bdraw[:,i + 1]);
        btt = copy(bt[:,i]);
        Vtt = copy(Symmetric(reshape(Vt[:,i], K, K)));
        fa = Symmetric(Vtt + Qt);
    
        invf = invss(fa);
        # invf = inv(f);
    
        cfe = bf - btt;
        bmean = btt + Vtt * invf * cfe;     # (Kx1)
        bvar = Symmetric(Vtt - Vtt * invf * Vtt);      # (KxK)
        bdraw[:,i] = randnx(bmean', bvar)';
    end

    return bdraw, bt, Vt
end
