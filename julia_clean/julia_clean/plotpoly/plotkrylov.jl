include("../geneigs.jl")
include("writefiles.jl")
include("../krylov.jl")

srand(5)

#eigdist(x)=exp(-(x-10.0)^2)+exp(-(x+10.0)^2)
eigdist(x)=exp(-(x-5)^2)
#eigdist(x)=1.0
#eigdist(x)=x+6
emin=-10.0
emax=10
n=100

gendeigs=geneigs(n,eigdist,emin,emax)

X=rand(n,n)
(Q,R)=qr(X)
X[:]=Q

shift=0.5*(gendeigs[convert(Int64,n/2)]+gendeigs[convert(Int64,n/2)+1])#4.6

trueeigs=(gendeigs-shift).^2
#trueeigs=gendeigs.^2
#trueeigs=1.0*gendeigs

A=X*diagm(trueeigs)*X'
A0=X*diagm(gendeigs)*X'

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
#=println("Est eigs:")
for i in 1:m*k
	println(xest[:,i]'*A0*xest[:,i])
end=#
Aq=xest'*(A0-shift*eye(n))*xest
Bq=xest'A*xest
#Aq=xest'*A0*xest
#Bq=xest'*xest
(lest2,xest2)=eig(Aq,Bq)
xest2=xest*xest2
println("A0 eigs:")
for i in 1:m*k
	#li=lest2[i]
	li=shift+(xest2[:,i]'*(A0-shift*eye(n))*xest[:,i])[1]
	
	println(li,"   ",norm(A0*xest2[:,i]-li*xest2[:,i])/norm(xest2[:,i]))
	#println(xest2[:,i]'*A0*xest2[:,i])
	#println(shift+(xest2[:,i]'*(A0-shift*eye(n))*xest[:,i])[1])
end

println("shift=$shift")

println("Real eigs:",gendeigs[1],"  ",gendeigs[n])
for i in 1:n
	#println(gendeigs[i])
end


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
