#using PyPlot

include("feast_lin.jl")
include("krylov.jl")
include("linsolvers.jl")
include("geneigs.jl")
#importall juliaFEAST

srand(3)

n=1000

distmin=-50.0
distmax=21.0

#eigdist(x)=exp(x)
#eigdist(x)=exp(-10*(x-1)^2)
eigdist(x)=1.0
trueeigs=geneigs(n,eigdist,distmin,distmax)

bdist(x)=1
beigs=geneigs(n,bdist,1,2)
bx=rand(n,n)
(bx,r)=qr(bx)
lambda=diagm(beigs)
lambda=sqrt(lambda)
F=bx*lambda

F=eye(n)
B=F*F'

lambda=diagm(trueeigs)
lambda=convert(Array{Float64,2},lambda)
xacc=rand(n,n)
(xacc,r)=qr(xacc)
A=xacc*lambda*xacc'

A=F*A*F'

m0=1
nc=4
eps=1e-16
maxit=300

#linear solver function:
function usegmres(A,b)
	P=eye(n)
	restarts=20
	ksize=1
	X=gmres(A,b,P,restarts,ksize)

	#println("    lin sys res = ",norm(b-A*X)/norm(b))
	
	return X
end

function useCG(A,b)
	its=20
	P=eye(n)
	A2=A'*A
	b2=A'*b
	#X=regularCG(A2,b2,its)
	#X=blockCG(A2,b2,P,its)
	#X=blockCGbasic(A,b,its)
	X=blockCGbasicnorm(A,b,its)

	#println("    lin sys res = ",norm(b-A*X)/norm(b))
	
	return X
end

function useBCG(A,b)
	its=5
        P=eye(n)
        X=blockBCG(A,b,P,its)

        #println("    lin sys res = ",norm(b-A*X)/norm(b))

        return X
end

#linsolve(A,b)=\(A,b)
linsolve(A,b)=usegmres(A,b)
#linsolve(A,B)=useCG(A,B)
#linsolve(A,B)=useBCG(A,B)

eigmin=501#47
eigmax=575#54

if eigmin>1
	emin=0.5*(trueeigs[eigmin]+trueeigs[eigmin-1])
else
	emin=0.5*(distmin+trueeigs[1])
end

if eigmax<n
	emax=0.5*(trueeigs[eigmax]+trueeigs[eigmax+1])
else
	emax=0.5*(distmax+trueeigs[eigmax])
end

x0=rand(Float64,n,m0)

emid=((emax+emin)/2.0)
demid=((emax-emin)/2.0)
emin2=0.0
emax2=demid^2
A2=(A-emid*eye(n))^2

A2=A
emin2=emin
emax2=emax



#(E,X,res,its)=feast_lin(A2,B,x0,nc,emin2,emax2,m0,eps,maxit,linsolve)

#(E,X)=basic_krylov(A2,x0,50)
(E,X)=folded_krylov(A2,x0,200,emin,emax)
#(E,X)=restart_krylov(A2,x0,20,20)
#(E,X)=mid_krylov(A,x0,100)

println("\n [Emin,Emax] = [ $emin, $emax]\n")

#for i in 1:its
#	println("     $i   $(res[i])")
#end

println("\n\nMeasured eigs:")
#println(real(sort(E)))
println("\n\nActual Eigs:")
#println(real(trueeigs[eigmin:eigmax]),"\n\n")

#println(diag(lambda))

#plot(res)

rx=A*X-B*X*diagm(E)
#rx=A*X-B*X*X'*A*X*inv(X'*B*X)
X2=A*X
println("Residuals:")

m02=size(X,2)
resE=zeros(m02,2)

for i in 1:m02
	resE[i,1]=norm(rx[:,i])/norm(X[:,i])
	resE[i,2]=E[i]
end

sortedResE=sortrows(resE,by=x->x[1])

for i in 1:m02
	println("$i    $(sortedResE[i,1])    $(sortedResE[i,2])")
end
