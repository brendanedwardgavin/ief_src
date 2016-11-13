#using PyPlot

#TRY USING FOLDED SPECTRUM

#include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
include("readMatrix.jl")
#importall juliaFEAST

srand(3)

#=
n=1000

distmin=-51.0
distmax=53.0

#eigdist(x)=exp(x)
eigdist(x)=exp(-10*(x-1)^2)
#eigdist(x)=1.0
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
xacc=rand(Complex128,n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'

A=F*A*F'

#A=A'*A
=#

############################################
println("Reading matrix.")
A=dreadMatrixCoo("../matrices/Na5_A.mtx")
#A=dreadMatrixCoo("../matrices/parsec/Si2-7e2/Si2.mtx")

n=size(A,1)

println("Reading done. n = $n.\n  Symmeterizing matrix.")

A=A+A'
A=A-0.5*spdiagm(diag(A))

println("Symmeterizing done. Doing calculation.")

##############################################

A=convert(SparseMatrixCSC{Complex128,Int64},A)

m=2

#b=rand(Float64,n,m)
#b=rand(Complex128,n,m)
#b=xacc[1:n,1:m]
bs=zeros(Complex128,n,m)
bs[1:4,1]=[1.0,2.0,3.0,4.0]
if m>1 
	bs[5:8,2]=[1.0,2.0,3.0,4.0]
end
b=A'*bs

#X=regularCG(A,b,10)
#X=blockCG(A,b,inv(A),5)
#X=blockBCG(A,b,eye(n),20)
#X=gmres(A,b,A,5,1)

X=gmres4(A,b,10,5,1e-16)
println("Sol=",X[1,1])

#X=blockCGbasicnorm(A,b,2500)
#X=blockCGbasic(A,b,50)
#X=stationaryIt(A,b,2,10)
#X=fom(A,b,40,1e-16)
#X1=affine(A,b,n*10)
#err1=norm(b-A*X1)/norm(b)
#println("err1=$err1\nerr2=$err2\nerr3=$err3\n")

println("\n\nLin sys residual = ",norm(b-A*X)/norm(b),"\n\n")
