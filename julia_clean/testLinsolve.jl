#using PyPlot

#TRY USING FOLDED SPECTRUM

#include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
#importall juliaFEAST

srand(3)

n=1000

distmin=51.0
distmax=53.0

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

m=1

b=rand(n,m)
#b=xacc[1:n,1:m]

#X=regularCG(A,b,10)
#X=blockCG(A,b,inv(A),5)
#X=blockBCG(A,b,eye(n),20)
#X=gmres(A,b,A,5,1)

#X=blockCGbasicnorm(A,b,2500)
#X=blockCGbasic(A,b,50)
#X=stationaryIt(A,b,2,10)
#X=fom(A,b,40,1e-16)
X1=affine(A,b,n*10)
err1=norm(b-A*X1)/norm(b)

X2=affine(A,b,n*10)
err2=norm(b-A*X2)/norm(b)

X3=0.5*(X1+X2)
err3=norm(b-A*X3)/norm(b)

println("err1=$err1\nerr2=$err2\nerr3=$err3\n")

#println("\n\nLin sys residual = ",norm(b-A*X)/norm(b),"\n\n")
