
include("../geneigs.jl")
include("../feast_lin.jl")

#generate initial eigenvalue distribution

n=545
lmin=-1
lmax=1

function seigdist(x)
	if (x<=1 && x>=-1)
		return 50/2
	elseif x>1
		return 495/abs(20.81-1.01)
	else
		return 0
	end
end


l=geneigs(n,seigdist,-1,20.81)


#make diagonal matrix of eigenvalues

A=diagm(l)

#use one FEAST multiplication rho*identity to get new eigenvalue distribution

newX=rhomult(A,eye(n),nc,lmin,lmax)

for i in 1:n
	println(newX[i,i])
end 

#save old and new distributions as DAT files


#save contour and contour points too?
