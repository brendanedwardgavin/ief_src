include("linsolvers.jl")

srand(34)

n=100
m=3
k=3

A=rand(n,n)
x0=rand(n,m)

#(V,H)=blockArnoldi(A,x0,k)


(V,H,Bsm)=blockArnoldiStart(A,x0,k)
for i in 2:k
	blockArnoldiPlus(A,V,H,i,k,m)
end


println("H before:\n\n",1*(H.!=0))

println("Bsm:\n\n",1*(Bsm.!=0))
#println("\n\nV'AV:\n\n",V'*A*V)

for i in 1:k
	vsm=H[(i-1)*m+1:i*m+m,(i-1)*m+1:i*m]
	M=blockGivens(vsm)
	Mb=eye((k+1)*m)
	Mb[(i-1)*m+1:i*m+m,(i-1)*m+1:i*m+m]=M

	println("\n\nMb:\n",1*(abs(Mb).>1e-14))

	H=Mb*H
	Bsm=Mb*Bsm
end

println("\n\nH after:\n\n",1*(abs(H).>1e-14))

println("Bsm:\n\n",1*(abs(Bsm).>1e-16))
