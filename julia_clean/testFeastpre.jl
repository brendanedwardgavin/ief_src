#using PyPlot

#TRY USING FOLDED SPECTRUM

include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
#importall juliaFEAST

srand(3)

n=1000

distmin=0.0
distmax=10.0

#eigdist(x)=exp(x)
eigdist(x)=exp(-10*(x-1)^2)
#eigdist(x)=1.0
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

m0=10
nc=4
eps=1e-12
maxit=300

#linear solver function:
function usegmres(A,b,y)
	P=eye(n)
	#P=(eye(n)+y*(A*y)')
	#P=A
	restarts=5
	ksize=1
	X= gmres(A,b,P,restarts,ksize)

	#println("    lin sys res = ",norm(b-A*X)/norm(b))
	
	return X
end

function useCG(A,b,y)
	its=20
	P=eye(n)
	#P=inv(A)
	X=blockCG(A,b,P,its)

	#println("    lin sys res = ",norm(b-A*X)/norm(b))
	
	return X
end

function useBCG(A,b,y)
	its=5
        P=inv(A)#eye(n)
        X=blockBCG(A,b,P,its)

        #println("    lin sys res = ",norm(b-A*X)/norm(b))

        return X
end

#linsolve(A,b,y)=\(A,b)
#linsolve(A,b,y)=usegmres(A,b,y)
#linsolve(A,B,y)=useCG(A,B,y)
linsolve(A,B,y)=useBCG(A,B,y)

#eigmin=507
#eigmax=514

eigmin=47
eigmax=54

#eigmin=n-7
#eigmax=n

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

#emin2=0.0
#emax2=(emax-(emin+emax)/2.0)^2
#A2=(A-((emin+emax)/2.0)*eye(n))^2

A2=A
emin2=emin
emax2=emax


(E,X,res,its)=feast_lin_pre(A2,B,x0,nc,emin2,emax2,m0,eps,maxit,linsolve)

println("\n [Emin,Emax] = [ $emin, $emax]\n")

#for i in 1:its
#	println("     $i   $(res[i])")
#end

println("\n\nMeasured eigs:")
println(real(sort(E)))
println("\n\nActual Eigs:")
println(real(trueeigs[eigmin:eigmax]),"\n\n")

#println(diag(lambda))

#plot(res)

#rx=A*X-B*X*diagm(E)
rx=A*X-B*X*X'*A*X*inv(X'*B*X)
X2=A*X
println("Residuals:")
for i in 1:m0
	println("$i   $(norm(rx[:,i])/norm(X[:,i]))")
end
