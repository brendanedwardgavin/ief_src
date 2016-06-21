#using PyPlot

include("feast.jl")
#include("psolve.jl")
include("readEigs.jl")

#importall juliaFEAST

srand(3)

n=100
n2=50
lambda=diagm(1:n)
lambda=convert(Array{Float64,2},lambda)
lambda[3,3]=3.9
#lambda=lambda-50.0*eye(n)
#lambda=lambda*(-1.0)
#for i in 1:n2
#	lambda[i,i]=1.0*i
#	lambda[i+n2,i+n2]=1.0*i 
#end
#lambda[1,1]=-50.0
xacc=rand(n,n)
#xacc=rand(Complex128,n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'
#A=A-50.5*eye(n)


m0=3
emin=0.0
emax=3.95

#A=readEigs("system2.eigs")
#=println("reading matrix")
A=zreadMat("system2.mat")
println("done. doing problem")
(n,n)=size(A)=#

#error("stopped")

#=m0=20
emin=-40.0
emax=-1e-1=#

#x0=rand(Complex128,n,m0)
x0=rand(Float64,n,m0)
#x0=convert(Array{Complex128,2},x0)
#x0=rand(Complex128,n,m0)
#A=convert(Array{Complex128,2},A)

eps=1e-12
maxit=100

#A=convert(Array{Complex128,2},A)

Asp=sparse(A)

#(E,X,res)=feast(Asp,x0,emin,emax,m0,eps,maxit)
(E,X,res,linits)=invIt(Asp,x0,emin,emax,m0,eps,maxit)
#(E,X,res,linits)=feast(Asp,x0,emin,emax,m0,3,eps,maxit)
#(E,X,res)=zfeast(Asp,x0,emin,emax,m0,eps,maxit)

#=for i in 1:size(res,1)
	println("     $i   $(res[i])")
end=#

#(E,X,res)=psolveInv(full(Asp),x0,eps,maxit)

println(real(E))

#println(diag(lambda))

#plot(res)
