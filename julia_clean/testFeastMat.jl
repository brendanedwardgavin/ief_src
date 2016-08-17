#using PyPlot

#TRY USING FOLDED SPECTRUM

include("feast_lin.jl")
include("linsolvers.jl")
include("geneigs.jl")
#importall juliaFEAST
include("readMatrix.jl")

println("Reading matrix.")
#A=dreadMatrixCoo("../matrices/Na5_A.mtx")
A=dreadMatrixCoo("../matrices/parsec/Si2-7e2/Si2.mtx")

n=size(A,1)

println("Reading done. n = $n.\n  Symmeterizing matrix.")

A=A+A'
A=A-0.5*spdiagm(diag(A))

println("Symmeterizing done. Doing calculation.")

n=size(A,1)

m0=60
nc=4
eps=1e-16
maxit=300
emin=17.052
emax=17.685

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
	its=10
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
	its=10
        P=eye(n)
        X=blockBCG(A,b,P,its)

        #println("    lin sys res = ",norm(b-A*X)/norm(b))

        return X
end

#linsolve(A,b)=\(A,b)
#linsolve(A,b)=usegmres(A,b)
linsolve(A,B)=useCG(A,B)
#linsolve(A,B)=useBCG(A,B)

x0=rand(Float64,n,m0)

emin2=0.0
emax2=(emax-(emin+emax)/2.0)^2
A2=(A-((emin+emax)/2.0)*eye(n))^2

A2=A
emin2=emin
emax2=emax

B=speye(n)

(E,X,res,its)=feast_lin(A2,B,x0,nc,emin2,emax2,m0,eps,maxit,linsolve)

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

#rx=A*X-B*X*diagm(E)
rx=A*X-B*X*X'*A*X*inv(X'*B*X)
X2=A*X
println("Residuals:")
for i in 1:m0
	#println("$i   $(norm(rx[:,i])/norm(X[:,i]))")
end
