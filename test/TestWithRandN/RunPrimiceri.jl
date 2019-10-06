
cd("/Users/zymantas/Dropbox/TVPVAR/_JuliaV1/TVPVARPkg/test/TestWithRandN")

using DelimitedFiles
using LinearAlgebra
include("extfunc.jl")
include("randnNEW.jl")
include("runtvp.jl")


# srand(5000)

data = readdlm("RunPrimiceri/primic_data.txt", ',');

mydraws = readdlm("mydraws.txt", ',');
global c = 0
global mydraws

REPS = 25
KEEP = 25
KEEPgap = 1

# REPS = 10000
# KEEP = 5000
# KEEPgap = 1

# REPS = 300000
# KEEP = 10000
# KEEPgap = 5

it_print = 1000
CONST = 1
P = 2
scaleb0var = 100
scaleQ0 = 0.0000000001

@time bmean, Qmean, sigmamean, bstore, Qstore, sigmastore, btstore, Vtstore = runtvpvarMLTVflat(data, REPS, KEEP, KEEPgap, it_print, CONST, P, scaleb0var, scaleQ0)

writedlm("RunPrimiceri/bmean.txt", bmean, ',');
writedlm("RunPrimiceri/Qmean.txt", Qmean, ',');
writedlm("RunPrimiceri/sigmamean.txt", sigmamean, ',');

    writedlm("StoredDraws/bstore2.txt", bstore, ',')
    writedlm("StoredDraws/Qstore2.txt", Qstore, ',')
    writedlm("StoredDraws/sigmastore2.txt", sigmastore, ',')
    writedlm("StoredDraws/btstore2.txt", btstore, ',')
    writedlm("StoredDraws/Vtstore2.txt", Vtstore, ',')
