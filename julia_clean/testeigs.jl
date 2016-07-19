include("geneigs.jl")

emin=-5.0
emax=5.0
neigs=1000

f(x)=exp(-10*(x-1)^2)+exp(-10*(x+1)^2)
#f(x)=exp(-10*(x-1)^2)
(a,err)=quadgk(f,emin,emax)
eigdist(x)=f(x)/a

eigs=geneigs(neigs,eigdist,emin,emax)

nsamples=neigs

#plt[:hist](eigs, bins=50)
