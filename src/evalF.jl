
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
function useN2fit(x::AbstractArray{Float64,3}, v::AbstractMatrix)
    q, w, e = size(x)
    y = similar(x, q, w)

    for a in 1:q, s in 1:w
        z = fit(Normal, x[a,s,:])
        y[a,s] = pdf.(z, v[a,s])
    end

    return y
end

# fitting the multivariate version
function useMvN2fit(x::AbstractArray{Float64,3}, v::AbstractMatrix)
    # w is no of variables
    # q is no of forecast horizon
    # e are number of draws
    q, w, e = size(x)
    y = similar(x, q)

    for a = 1:q
        z = fit(MvNormal, x[a,:,:])
        y[a] = pdf(z, v[a,:])
    end

    return y
end

#############################################
## fitting PIT 
"
PSPIT Calculate Probability Integral Transform time series.
   [Z, ZHIST] = PSPIT(FORECAST, OBS) obtains the PIT
   z-series for a given vector of observations in OBS, with associated
   empirical probability forecasts given in the matrix FORECASTS. Each
   row of FORECASTS contains an empirical distribution for each element
   of OBS. [...] = PSPIT(..., ZBINS) specifies the number of histogram bins
   on the interval [0,1] (default is 10 bins). [...] =
   PSPIT(..., ZBINS, ZSIGNIF) specifies a significance probability for the
   null hypothesis test of the z-series bin counts (default 0.95 - 95%
   significance level).
   Output parameter Z is the z-series, and ZHIST is a structure
   with the following elements:
    'h.counts'   - Count of the number of z-series elements falling within
                 each histogram bin on the interval [0,1].
    'h.centres'  - Z-series histogram bin centres.
    'confmean' - Expected z-series histogram count, under the null
                 hypothesis of z-series being iid uniform on [0,1]. 
    'confpos'  - Upper confidence interval for the z-series histogram
                 counts, under the null hypothesis of z-series being iid
                 uniform on [0,1].
    'confneg'  - Lower confidence interval for the z-series histogram
                 counts, under the null hypothesis of z-series being iid
                 uniform on [0,1].
(cc) Max Little, 2008. This software is licensed under the
Attribution-Share Alike 2.5 Generic Creative Commons license:
http://creativecommons.org/licenses/by-sa/2.5/
If you use this work, please cite:
Little MA et al. (2008), Parsimonious Modeling of UK Daily Rainfall for
Density Forecasting, in Geophysical Research Abstracts, EGU General
Assembly, Vienna 2008, Volume 10.
"
function pspit(forecasts::Matrix, obs::Vector; zsignif::Float64 = 0.95, nbins::Int = 10, returnZ::Bool=false)
    N = length(obs)
    M, S = size(forecasts)

    @assert M == N "Length of OBS must match number of rows in FORECASTS"

    z = similar(obs)
    freqs = (1. : S) ./ S
    scale = similar(forecasts, S)
    for i = 1:N
        scale .= sort(forecasts[i,:])

        @views obsI = obs[i]
        if obsI < scale[1]
            z[i] = 0.
        elseif obsI > scale[end]
            z[i] = 1.
        else
            idx = findlast(scale .<= obsI)
            if idx > 1
                z[i] = freqs[idx - 1] + rand(Uniform()) * (freqs[idx] - freqs[idx - 1])
            else
                z[i] = 0.
            end
        end
    end

    if returnZ
        return z
    end

    # Obtaining a histogram
    bins = LinRange(0., 1., nbins+1)
    h = fit(Histogram, z, bins)

    # Find the approximate confintv confidence interval for
    # binomially-distributed PIT histogram bin counts, under
    # the assumption of iid U(0,1) PIT z-series

    bprob = 1. / nbins
    meanheight = convert(Int, floor(N/nbins))::Int
    binocumul(x) = cdf(Binomial(N,bprob), x)
    for confwidth = 1:meanheight
        cip = meanheight + confwidth
        cin = meanheight - confwidth
        ci = binocumul(cip) - binocumul(cin - 1.)
        if ci >= zsignif
            return (z=z, zhist=h, confmean= meanheight, confpos = cip, confneg = cin)
        end
    end

    @warn "the confidence interval could not be compute due to few observation"
    return (z=z, zhist=h, confmean= meanheight)

end