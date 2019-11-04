module myJuliaUtils

using Random
using LinearAlgebra
using Statistics
using Distributions

export

# extfunc.jl
lag0, eye, vech, vec2sym, vec2ltri, vec2ltriW1, stabcheck, stabcheckC, preparexy, sumsqr, nanmean,
wish, iwish, cholPSD, randnPSD, genPSDmatStrict, genPDmat, genPSDmat, simMC, quantileArr, companionf, 
ismyapprox, normpdf, acf, decVCV, ols1, VARols1, inbetween,
# evalF.jl
useN2fit, useMvN2fit,
# data_manipulation.jl
dataTransf, detrend,
# generate.jl
genLPM

include("extfunc.jl");
include("evalF.jl");
include("data_manipulation.jl");
include("generate.jl");



end # module
