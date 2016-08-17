using PyPlot
include("geneigs.jl")

emin=-50.0
emax=10.0
neigs=1000

#f(x)=exp(-10*(x-1)^2)+exp(-10*(x+1)^2)
#f(x)=exp(-10*(x-1)^2)
f(x)=1.0

eigs=geneigs(neigs,f,emin,emax)

nsamples=neigs

plt[:hist](eigs, bins=50)
