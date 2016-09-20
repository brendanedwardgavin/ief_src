#using PyPlot

include("feast.jl")
include("readEigs.jl")
#importall juliaFEAST

srand(3)

n=100
lambda=diagm(1:n)
lambda=convert(Array{Float64,2},lambda)
#lambda[3,3]=3.9

xacc=rand(n,n)
#xacc=rand(Complex128,n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'
#A=A-50.5*eye(n)


m0=10
emin=50.5
emax=53.95
mactual=3

#A=readEigs("system2.eigs")
#=println("reading matrix")
A=zreadMat("system2.mat")
println("done. doing problem")
(n,n)=size(A)=#

x0=rand(Float64,n,m0)

eps=1e-12
maxit=200

#A=convert(Array{Complex128,2},A)

(E,X,res,its)=feast(A,x0,emin,emax,m0,mactual,eps,maxit)

for i in 1:its
	println("     $i   $(res[i])")
end

println(real(E))

#println(diag(lambda))

#plot(res)
