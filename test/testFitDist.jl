# using KernelDensity
using Distributions
using myJuliaUtils
using Test

# univariate case

# reps=200
reps=20000000
x = randn(reps)
a=0:1
# @time y = kde(x)
# pdf(y, a)

@time z=fit(Normal, x)
a1=pdf.(z,a)

@test isapprox(a1,pdf.(Normal(0.,1.),a); rtol=1e-3)



q=3
w=4
e=reps
x=randn(q,w,e);
v=reshape(collect(1:1:q*w) ./10,q,w)
o= sum(v.-0.1 .<=x .<=v.+0.1,dims=3)./e ./0.2
# @time y=TVPVARPkg.useKDE2fit(x,v)
# @btime y=findKD($x,$v)

@time z=useN2fit(x,v)
@test isapprox(z, pdf.(Normal(0.,1.),v) ; rtol=1e-3)


# multivariate case
@time z=useMvN2fit(x,v) 

f=similar(x,q)
for a=1:q
    f[a]=pdf(MvNormal(zeros(w),eye(w)), v[a,:])
end

@test isapprox(z, f ; rtol=1e-3)


