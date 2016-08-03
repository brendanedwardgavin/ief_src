#using PyPlot

#TRY USING FOLDED SPECTRUM

#include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
include("readMatrix.jl")
#importall juliaFEAST

println("Reading matrix.")
#A=dreadMatrixCoo("../matrices/Na5_A.mtx")
A=dreadMatrixCoo("../matrices/smallTest_A.mtx")

n=size(A,1)

println("Reading done. n = $n.\n  Symmeterizing matrix.")

A=A+A'
A=A-0.5*spdiagm(diag(A))

println("Symmeterizing done. Doing calculation.")

m=1

#b=rand(n,m)
b=zeros(n,m)
b[1,1]=1.0
b[2,1]=2.5

#X=regularCG(A,b,10)
#X=blockCG(A,b,inv(A),5)
#X=blockBCG(A,b,eye(n),20)
#X=gmres(A,b,A,5,1)

X=blockCGbasicnorm(A,b,2)
#X=blockCGbasic(A,b,10)

println("\n\nLin sys residual = ",norm(b-A*X)/norm(b),"\n\n")


println(A)
