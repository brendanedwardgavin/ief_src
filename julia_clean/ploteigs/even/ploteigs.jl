
include("../geneigs.jl")
include("../feast_lin.jl")

#generate initial eigenvalue distribution

#n=545
n=10
nin=4
eigmax=5
lmin=0
density=n/(eigmax-lmin)
lmax=lmin+nin/density
nc=8

println(lmax," ",eigmax)

#eigmax=1.1
#eigmax=20.81



#nin=50

function seigdist(x)

	if (x<=lmax && x>=lmin)
		#return 50/2
		return density
	elseif x>lmax
		#return 495/abs(eigmax-1.01)
		return density
	else
		return 0
	end
end

l=geneigs(n,seigdist,lmin,eigmax)
#lin=geneigs(nin,seigdist,lmin,lmax)
#lin=geneigs(nin,seigdist,-1,1)
#lout=geneigs(n-nin,seigdist,lmax,eigmax)

#l=[lin;lout]

nin=0

for i in 1:n
	if l[i]>=lmin && l[i]<=lmax
		nin=nin+1
	end
end

println("nin=$nin")

#make diagonal matrix of eigenvalues

A=diagm(l)

#use one FEAST multiplication rho*identity to get new eigenvalue distribution

newX=rhomult(A,eye(n),nc,lmin,lmax)

newX=-1.0*newX

#save old and new distributions as DAT files

open("xin.dat","w") do file
	for i in 1:nin
		write(file,"$(l[i])  0.0\n")
	end
end

open("xout.dat","w") do file
	for i in nin+1:n
		write(file,"$(l[i])  0.0\n")
	end
end


open("newxin.dat","w") do file
	for i in 1:nin
		write(file,"$(newX[i,i])  0.0\n")
	end
end


open("newxout.dat","w") do file
	for i in nin+1:n
		write(file,"$(newX[i,i])  0.0\n")
	end
end

#save contour and contour points too?
