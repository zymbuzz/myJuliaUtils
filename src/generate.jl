function genLPM(lammbda::AbstractVector, T::Int, X::AbstractMatrix)
    # generating the Linear Probability Model
    # lammbda is a vector (n x 1)
    # X is the matrix (T x n)

    y = BitArray(undef, T)
    for i = 1:T
        y[i] = dot(X[i,:], lammbda) > rand()
    end

    return y

end

genMC(PIt::Array{Float64,3}, T::Int; Tdrop::Int = 0, S0::Int = 1) = genEndoMC(PIt, T; Tdrop = Tdrop, S0 = S0)

function genEndoMC(PIt::Array{Float64,3}, T::Int; Tdrop::Int = 0, S0::Int = 1)
    S  = [S0]; # fix starting point
    for i = 2:T + Tdrop
        cumu_p = cumsum(PIt[S[i - 1], :,i])
        push!(S, sum(rand() .> cumu_p) + 1)
    end
    return S[Tdrop + 1:end]
end

simMC(PI::Array{Float64,2}, T::Int; Tdrop::Int = 0, S0::Int = 1) = genMC(PI, T, Tdrop = Tdrop, S0 = S0)

function genMC(PI::Array{Float64,2}, T::Int; Tdrop::Int = 0, S0::Int = 1)
    S  = [S0] # fix starting point
    cumPI = cumsum(PI;dims = 2)
    for i = 2:T
        push!(S, sum(rand() .> cumPI[S[i - 1], :]) + 1);
    end
    return S[Tdrop + 1:end]
end