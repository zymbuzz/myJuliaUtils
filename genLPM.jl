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
