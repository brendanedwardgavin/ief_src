include("linsolvers.jl")

srand(34)

n=100
m=2
k=3

A=rand(n,n)
x0=rand(n,m)

#(V,H)=blockArnoldi(A,x0,k)


(V,H,Bsm)=blockArnoldiStart(A,x0,k)
for i in 2:k
	blockArnoldiPlus(A,V,H,i,k,m)
end


println("H:\n\n",1*(H.!=0))

#println("Bsm:\n\n",1*(Bsm.!=0))
#println("\n\nV'AV:\n\n",V'*A*V)
