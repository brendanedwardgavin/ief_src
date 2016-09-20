#using PyPlot

include("feast.jl")
include("readMatrix.jl")

A=dreadMatrixCoo("test.mtx")

(n,m)=size(A)

m0=3
emin=3.5
emax=5.5
mactual=2

nc=16

#A=readEigs("system2.eigs")
#=println("reading matrix")
A=zreadMat("system2.mat")
println("done. doing problem")
(n,n)=size(A)=#

x0=rand(Float64,n,m0)

eps=1e-12
maxit=200

#A=convert(Array{Complex128,2},A)

(E,X,res,its)=feast(A,x0,nc,emin,emax,m0,eps,maxit)

#for i in 1:its
#	println("     $i   $(res[i])")
#end

println(real(E))

#println(diag(lambda))

#plot(res)
