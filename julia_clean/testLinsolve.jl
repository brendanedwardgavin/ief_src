#using PyPlot

#TRY USING FOLDED SPECTRUM

#include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
#importall juliaFEAST

srand(3)

n=1000

distmin=-50.0
distmax=3.0

#eigdist(x)=exp(x)
#eigdist(x)=exp(-10*(x-1)^2)
eigdist(x)=1.0
trueeigs=geneigs(n,eigdist,distmin,distmax)

bdist(x)=1
beigs=geneigs(n,bdist,1,2)
bx=rand(n,n)
(bx,r)=qr(bx)
lambda=diagm(beigs)
lambda=sqrt(lambda)
F=bx*lambda

F=eye(n)
B=F*F'

lambda=diagm(trueeigs)
lambda=convert(Array{Float64,2},lambda)
xacc=rand(n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'

A=F*A*F'

#A=A'*A

m=5

b=rand(n,m)

#X=regularCG(A,b,10)
#X=blockCG(A,b,inv(A),5)
#X=blockBCG(A,b,eye(n),20)
#X=gmres(A,b,A,5,1)

#X=blockCGbasicnorm(A,b,2500)
X=blockCGbasic(A,b,2500)


println("\n\nLin sys residual = ",norm(b-A*X)/norm(b),"\n\n")
