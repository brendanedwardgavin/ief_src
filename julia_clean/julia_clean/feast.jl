

#module juliaFEAST

#export zfeast,feast,feastG,invIt,alinvIt

#include("gmresRP.jl")
include("lsSolve.jl")
#include("psolve.jl")

using FastGaussQuadrature

gauss_weight8=[0.362683783378361;0.362683783378361;0.313706645877887;0.313706645877887;0.222381034453374;0.222381034453374;0.101228536290376;0.101228536290376]
gauss_ab8=[0.183434642495649;-0.183434642495649;0.525532409916328;-0.525532409916328;0.796666477413626;-0.796666477413626;0.960289856497536;-0.960289856497536]

function feastfilt(A,Y,lmin,lmax)
	nc=8
	(n,m)=size(Y)

	r=(lmax-lmin)/2.0
	
	
	Q=zeros(n,m)	
	
	for k in 1:nc
		theta=-1.0*(pi/2)*(gauss_ab8[k]-1)
		zk=(lmin+lmax)/2+r*exp(0+im*theta)
		G=zk*speye(n)-A
		Qk=\(G,Y)
		Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
	end

	return Q	
end

function invIt(A::SparseMatrixCSC{Float64,Int64},x0::Array{Float64,2},lmin::Float64,lmax::Float64,m0::Int64,eps::Float64,maxit::Int64)
	
	n::Int64=size(x0,1)
	res::Float64=1.0

	x=x0;

	lest=zeros(Float64,m0,1)
	lest=diag(x0'*(A*x0))

	it=0

	maxlinit=3000

	residuals=zeros(Float64,maxit)

	lest=zeros(m0)
	reslist=zeros(m0)
	linits=zeros(maxit)
	while res>eps && it<maxit
		it=it+1
		Q=zeros(x)	
		#println(size(Q))
				
		restarts=5000
		subspaces=10

		Gr=convert(Array{Complex128,2},full(A)-lmin*eye(n))
		br=convert(Array{Complex128,2},x)
		#Q=\(A-lmin*eye(n),x)

		Pi=zeros(n,n)	
		#P=diagm(diag(Gr))
		#P=eye(n)
		srand(1)
		randmat=rand(n,n)
		P=Gr+3.0*norm(Gr)*randmat/norm(randmat)
		println("   ILU error=$(norm(Gr-P)/norm(Gr))")
		Pi=Gr*br*br'+P*(eye(n)-br*br')		
		#Pi=br*br'*Gr'+eye(n)
		#Pi=eye(n)
		#Pi=P
		
		#Pi=eye(n)#br*br'*Gr'+eye(n)
		(Q,linits[it],linres)=improveSolvePeps(Gr,br,restarts,subspaces,1e-1*min(res,1),Pi)
		#(Q,linits[it],linres)=improveSolveSub(Gr,br,restarts,subspaces,1e-3*min(res,1),Pi)

		#error("done")
		#Q=improveSolvePeps(Gr,br,restarts,subspaces,1e-2)
		
		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*Q)		
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		
		resvecs=A*x-x*spdiagm(lest)
		reslist=zeros(m0,1)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		(res,ind)=findmax(reslist)
		
		println("invIt Iteration $it: res = $res")
		#println(lest[1:5])
		#println("~~~X norm = $(norm(x))")
		if isnan(res)
			#println("Got NaN! Printin X...")
			#println(x)
			break
		end
		residuals[it]=res
		#println(lest)
		#println(res)
	end

	println("Final residuals: $reslist")
	return (lest,x,residuals,linits)
end

function zfeast(A::SparseMatrixCSC{Complex128,Int64},x0::Array{Complex128,2},lmin::Float64,lmax::Float64,m0::Int64,eps::Float64,maxit::Int64)
	nc::Int64=8
	n::Int64=size(x0,1)
	res::Float64=1.0

	r::Float64=(lmax-lmin)/2.0
	x=x0;

	lest=zeros(Float64,m0,1)
	lest=diag(x0'*(A*x0))

	it=0


	maxlinit=3000

	residuals=zeros(Float64,maxit)

	lest=zeros(m0)
	while res>eps && it<maxit
		it=it+1
		Q=zeros(x)	
		#println(size(Q))
		
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_ab8[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			G=zk*speye(n)-A	
			
			Qk=zeros(Complex128,n,m0)
			Gr=full(G)
			br=convert(Array{Complex128,2},x)
			
			restarts=5
			subspaces=2
			#lest=[98.0,99.0,100.0]			
			#=if it>1
				#x0=x
				#x0=x*inv(diagm(lest))
				x0=zeros(n,m0)
				Qk=improveSolveIn(Gr,br,x0,restarts,subspaces)
			else
				Qk=improveSolve(Gr,br,restarts,subspaces)
			end=#
	
			x0=zeros(n,m0)
			#Qk1=improveSolveIn(Gr,br,x0,restarts,subspaces)	

			x0=zeros(n,m0)
			#Qk2=improveSolveIn(Gr',br,x0,restarts,subspaces)

			Qk1=\(G,x)
			Qk2=\(G',x)			

			p=1	
			delta=norm(G*Qk[:,p]-x[:,p])/norm(x[:,p])
			#println("        Lin sys error = $delta")
			
			#Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
		
			Q=Q-0.25*gauss_weight8[k]*r*(exp(im*theta)*Qk1+exp(-1.0*im*theta)*Qk2)
		end

		#(QQ,RR)=qr(Q)
		#Q[:,1:m0]=QQ[:,1:m0]
		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*Q)		
		(u,s,v)=svd(Bq)
		#println(s)
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		
		resvecs=A*x-x*spdiagm(lest)
		reslist=zeros(m0,1)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		(res,ind)=findmin(reslist)
		
		println("FEAST Iteration $it: res = $res")
		#println(real(lest))
		#println("~~~X norm = $(norm(x))")
		if isnan(res)
			#println("Got NaN! Printin X...")
			#println(x)
			break
		end
		residuals[it]=res
		#println(lest)
		#println(res)
	end
	return (lest,x,residuals)
end

function feast(A,x0,nc,lmin,lmax,m0,eps,maxit)
	#nc::Int64=8
	n::Int64=size(x0,1)
	res::Float64=1.0

	r::Float64=(lmax-lmin)/2.0
	x=x0;

	lest=zeros(Float64,m0,1)
	lest=diag(x0'*(A*x0))

	it=0

	residuals=zeros(Float64,maxit)

	lest=zeros(m0)
	linits=zeros(maxit)
	linres=zeros(maxit)
	
	resvecs=zeros(Float64,n,m0)

	(gauss_point,gauss_weight)=gausslegendre(nc)
	while res>eps && it<maxit
		it=it+1
		Q=zeros(x)	
		#println(size(Q))
		linresnc=zeros(nc)
		linitsnc=zeros(nc)	
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_point[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			G=zk*speye(n)-A
			
			Qk=\(G,x)
			
			restarts=5
			krylov=1
	
			#if it>1 
			#	Qk=improveSolveR(G,x,restarts,krylov,resvecs)
			#else
			#	Qk=improveSolve(G,x,restarts,krylov)
			#end
			#Qk=improveSolve(G,x,restarts,krylov)
			Q=Q-(gauss_weight[k]/2)*real(r*exp(im*theta)*Qk)
		end

		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*Q)		
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		F=eigfact(Bq)
		best=F[:values]
		mactual=0
		for i in 1:m0
			if best[i]>= 0.0 && lest[i] >= lmin && lest[i] <=lmax
				mactual=mactual+1
			end
		end
		
		resvecs=A*x-x*spdiagm(lest)
		reslist=zeros(m0)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		#(res,ind)=findmax(reslist)
		reslistsort=sort(reslist)
		res=reslistsort[mactual]		
			
		println("FEAST Iteration $it: res = $res,  m=$mactual")
		#println("    Lin sys res=$(linres[it]), its=$(linits[it])")
		
		if isnan(res)
			println("Got NaN! Printin X...")
			#println(x)
			break
		end
		residuals[it]=res
		#println(lest)
		#println(res)
	end
	return (lest,x,residuals,it)
end

function feastG(A::SparseMatrixCSC{Float64,Int64},B::SparseMatrixCSC{Float64,Int64},x0::Array{Float64,2},lmin::Float64,lmax::Float64,m0::Int64,eps::Float64,maxit::Int64)
	nc::Int64=8
	n::Int64=size(x0,1)
	res::Float64=1.0

	r::Float64=(lmax-lmin)/2.0
	x=x0;

	lest=zeros(Float64,m0,1)
	lest=diag(x0'*(A*x0))

	it=0


	maxlinit=3000

	residuals=zeros(Float64,maxit)

	while res>eps && it<maxit
		it=it+1
		Q=zeros(x)	
		#println(size(Q))
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_ab8[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			G=zk*B-A

			Bf=full(B)
			Af=full(A)
			Gf=zk*Bf-Af			

			#Qk=\(G,B*x)
			#Gf=full(G)
			xr=convert(Array{Complex128,2},x)
			
			#Qk=improveSolveG(Gf,Bf,xr,5,1)

			restarts=5
                               
			#x0=zeros(n,m0)
                        #Qk=improveSolveGIn(Gf,Bf,xr,x0,restarts,1)
			Qk=\(Gf,Bf*xr)                      

			p=1	
			delta=norm(G*Qk[:,p]-B*x[:,p])/norm(B*x[:,p])
			#println("        Lin sys error = $delta")
			
			Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
		end

		#(QQ,RR)=qr(Q)
		#Q[:,1:m0]=QQ[:,1:m0]
		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*(B*Q))			
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		
		resvecs=A*x-B*x*spdiagm(lest)
		reslist=zeros(m0,1)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(B*x[:,i])
		end
		(res,ind)=findmin(reslist)
		
		println("FEAST Iteration $it: res = $res")
		#println("~~~X norm = $(norm(x))")
		if isnan(res)
			#println("Got NaN! Printin X...")
			#println(x)
			break
		end
		residuals[it]=res
		#println(lest)
		#println(res)
	end
	return (lest,x,residuals)
end

#end
