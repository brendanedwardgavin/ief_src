#using PyPlot
include("lsSolve.jl")
include("gmresCRP.jl")

srand(1)

n=100
m=2

#A=2*rand(n,n)-ones(n,n)
#A=A+A'

lambda=diagm(1:n)
#lambda=convert(Array{Float64,2},lambda)
lambda=convert(Array{Complex128,2},lambda)
for i in 1:n
	lambda[i,i]=lambda[i,i]+5*im
end
xacc=rand(n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'

B=2*rand(n,m)-ones(n,m)
B=convert(Array{Complex128,2},B)
errs=zeros(n)

X=zeros(Complex128,n,m)

#X=improveSolve(A,B,1,10)
#X=lsSolve(A,B,20)

x0=zeros(n,1)
m0=10
for i in 1:m
	B1=zeros(Complex128,n,1)
	B1[:,1]=B[:,i]
	X[:,i]=gmresCRP(A,x0,B[:,i],eye(n),1,m0,1e-17)	
	#X[:,i]=lsSolve(A,B1,m0)
	#X[:,i]=improveSolve(A,B1,1,m0)	
end

println("error=",norm(B-A*X)/norm(B))
#=
for m0 in 1:20#n-1
	println(m0) 
	#X=lsSolve(A,B,m0)

	x0=zeros(n,1)
	for i in 1:m
		B1=zeros(n,1)
		B1[:,1]=B[:,i]
		X[:,i]=gmresRP(A,x0,B[:,i],eye(n),10,m0,1e-17)	
		#X[:,i]=lsSolve(A,B1,m0)
		#X[:,i]=improveSolve(A,B1,2,m0)	
	end
	#X=improveSolve(A,B,2,m0)
	errs[m0]=norm(B-A*X)/norm(B)
end

plot(errs)=#
