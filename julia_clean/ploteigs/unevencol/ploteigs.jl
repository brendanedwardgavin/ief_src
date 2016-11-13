
include("../../geneigs.jl")
include("../../feast_lin.jl")

#generate initial eigenvalue distribution

#n=545
n=26
nin=6

eigmax=5
eigmin=0
density=10/(eigmax-eigmin)

lmax=eigmax
lmin=4/density
#lmax=lmin+nin/density
nc=3

inmin=n-nin+1
inmax=n-nin+2

println(lmax," ",eigmax)

#eigmax=1.1
#eigmax=20.81
#nin=50

function seigdist(x)
	return density
end

function deigdist(x)
	return density*x^2
end

#l=geneigs(n,seigdist,eigmin,eigmax)
#lin=geneigs(nin,seigdist,lmin,lmax)
lin=geneigs(10,seigdist,eigmin,eigmax)
lminfeast=(lin[5]+lin[4])/2
lmaxfeast=(lin[6]+lin[7])/2

#lout=geneigs(n-nin,deigdist,eigmin,lminfeast)#lmin)
lout=geneigs(n-nin,deigdist,eigmin,lin[5])
l=[lout;lin[5:10]]

dl=lout[n-nin]-lin[5]
lminfeast=(lout[n-nin]+lin[5])/2+0.45*dl
 

nin=0

for i in 1:n
	if l[i]>=lminfeast && l[i]<=lmaxfeast
		nin=nin+1
	end
end

println("nin=$nin")

#make diagonal matrix of eigenvalues

A=diagm(l)

#use one FEAST multiplication rho*identity to get new eigenvalue distribution
#lminfeast=(l[inmin]+l[inmin-1])/2
#lmaxfeast=l[26]+1.1*l[26]#(l[inmax]+l[inmax+1])/2
newX=rhomult(A,eye(n),nc,lminfeast,lmaxfeast)

newX=-1.0*newX

#save old and new distributions as DAT files

open("xin.dat","w") do file
	for i in 1:n
		if (i>=inmin && i<=inmax)
		write(file,"$(l[i])  0.0\n")
		end
	end
end

open("xout.dat","w") do file
	for i in 1:n
		if !(i>=inmin && i<=inmax)
		write(file,"$(l[i])  0.0\n")
		end
	end	
end


open("newxin.dat","w") do file
	for i in 1:n
		if (i>=inmin && i<=inmax)
			write(file,"$(newX[i,i])  0.0\n")
		end
	end
end


open("newxout.dat","w") do file
	for i in 1:n
		if !(i>=inmin && i<=inmax)
			write(file,"$(newX[i,i])  0.0\n")
		end
	end
	
end

#save contour and contour points too?
