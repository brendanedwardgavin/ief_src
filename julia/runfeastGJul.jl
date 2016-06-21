#using PyPlot

include("feast.jl")
#include("lsSolve.jl")

#importall juliaFEAST

srand(2)

n=100


bx=rand(n,n)
(bx,r)=qr(bx)
lambda=diagm(rand(n))
lambda=sqrt(lambda)
F=bx*lambda
B=F*F'
#F=eye(n)
#B=eye(n)

lambda=diagm(1:n)
lambda=convert(Array{Float64,2},lambda)
lambda[3,3]=3.9
xacc=rand(n,n)
(xacc,r)=qr(xacc)

A=xacc*lambda*xacc'
A=F*A*F'

m0=3
x0=rand(n,m0)

emin=0.0
emax=2.0

eps=1e-12
maxit=400

#A=inv(B)*A
#B=eye(n)
Asp=sparse(A)
Bsp=sparse(B)
#(E,X)=eig(A,B)

#(E,X,res)=feastG(Asp,Bsp,x0,emin,emax,m0,eps,maxit)
(E,X,res,linits)=invItG(Asp,B,x0,emin,emax,m0,eps,maxit)
#for i in 1:size(res,1)
#	println("     $i   $(res[i])")
#end

#esort=sort(E)
#println(esort[1:3])

println(E)

#plot(res)

