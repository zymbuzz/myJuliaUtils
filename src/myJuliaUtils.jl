module myJuliaUtils

using Random
using LinearAlgebra
using Statistics
using Distributions

export

# extfunc.jl
lag0, eye, vech, vec2sym, vec2ltri, vec2ltriW1, stabcheck, stabcheckC, preparexy, nanmean,
wish, iwish, cholPSD, randnPSD, genPSDmatStrict, genPDmat, genPSDmat, simMC, quantileArr, companionf, 
ismyapprox, normpdf, acf, decVCV, ols1, inbetween,
# evalF.jl
useN2fit, useMvN2fit

include("extfunc.jl");
include("evalF.jl");


end # module
