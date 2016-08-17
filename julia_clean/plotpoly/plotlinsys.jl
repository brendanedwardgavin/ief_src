include("../geneigs.jl")
include("writefiles.jl")
include("../krylov.jl")

#eigdist(x)=exp(-(x-10.0)^2)+exp(-(x+10.0)^2)
eigdist(x)=1.0
#eigdist(x)=x+6
emin=-5.0
emax=4.0
n=20

trueeigs=geneigs(n,eigdist,emin,emax)

X=rand(n,n)
(Q,R)=qr(X)
X[:]=Q

A=X*diagm(trueeigs)*X'

m=1
k=3

x0=rand(n,m)

(lest,xest)=basic_krylov(A,x0,k)

#=
for i in 1:n
	println("$i $(trueeigs[i])")
end
=#

println(lest)

function truepoly(x)
	fx=1.0
	for i in 1:n
		fx=fx*(x-trueeigs[i])	
	end

	return fx
end

function kpoly(x)
	fx=1.0
	for i in 1:k
                fx=fx*(x-lest[i])
        end

        return fx
end


writefunction(truepoly,500,minimum(trueeigs),maximum(trueeigs),"truepolyf.dat")

writefunction(kpoly,500,minimum(trueeigs),maximum(trueeigs),"kpolyf.dat")

writezeros(trueeigs,"trueeigs.dat")

writezeros(lest,"approxeigs.dat")
