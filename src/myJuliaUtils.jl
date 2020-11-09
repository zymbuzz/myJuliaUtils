module myJuliaUtils

using Random
using LinearAlgebra
using Statistics
using Distributions
using StatsBase: counts, Histogram, fit

export

# extfunc.jl
lag0, eye, vech, vechRshp, vec2sym, vec2ltri, vec2ltriW1, stabcheck, stabcheckC, preparexy, sumsqr, nanmean,
wish, iwish, cholPSD, regMat2PD!, randnPSD, genPSDmatStrict, genPDmat, genPSDmat, quantileArr, companionf, 
ismyapprox, normpdf, acf, decVCV, ols1, VARols1, inbetween, getmultdiag!,transf1To,transf1Back,split_data,standardize,unstandardize, showprog,
# evalF.jl
useN2fit, useMvN2fit, pspit, pscrps,
# data_manipulation.jl
dataTransf, detrend,
# generate.jl
genLPM, genEndoMC, simMC, genMC

include("extfunc.jl");
include("evalF.jl");
include("data_manipulation.jl");
include("generate.jl");

end # module
