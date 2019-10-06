
# # This  function fits Kernel Density and evaluates density at v (Needs KernelDensity.jl)
# function useKDE2fit(x::AbstractArray{Float64,3},v::AbstractMatrix)
#     q,w,e=size(x)
#     y = similar(x,q,w)

#     for a=1:q, s=1:w
#         z = kde( x[a,s,:] )
#         y[a,s] = pdf(z, v[a,s])
#     end

#     return y
# end

# This  function fits Normal pdf and evaluates density at v (Needs Distributions.jl)
function useN2fit(x::AbstractArray{Float64,3},v::AbstractMatrix)
    q,w,e=size(x)
    y = similar(x,q,w)

    for a=1:q, s=1:w
        z = fit(Normal, x[a,s,:])
        y[a,s]=pdf.(z, v[a,s])
    end

    return y
end

# fitting the multivariate version
function useMvN2fit(x::AbstractArray{Float64,3},v::AbstractMatrix)
    # w is no of variables
    # q is no of forecast horizon
    # e are number of draws
    q,w,e = size(x)
    y = similar(x, q)

    for a=1:q
        z = fit(MvNormal, x[a,:,:])
        y[a]=pdf(z, v[a,:])
    end

    return y
end
