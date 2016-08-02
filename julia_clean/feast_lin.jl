include("linsolvers.jl")

using FastGaussQuadrature

function feast_lin(A,B,x0,nc,lmin,lmax,m0,eps,maxit,linsolve)
	
	#linsolve(A,b) ->returns solution to A*x=b

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

		linresnc=zeros(nc)
		linitsnc=zeros(nc)	
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_point[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			G=zk*B-A
			
			Qk=linsolve(G,B*x)#\(G,x)

			Q=Q-(gauss_weight[k]/2)*real(r*exp(im*theta)*Qk)
		end

		#=for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end
		=#

		#(Q,R)=qr(Q)

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*B*Q)		
		
		try
			F=eigfact(Hermitian(Aq),Hermitian(Bq))
		catch theerror
			F=eigfact(Bq)
			lest=F[:values]
			println("Bq eigenvalues:\n $(lest)\n")
			error("Error in reduced eigenvalue problem\n Original error:\n $theerror")
		end

		#(val,vec)=eig(Hermitian(Bq))
		#println("   $val")

		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		F=eigfact(Bq)
		best=F[:values]
		mactual=0
		for i in 1:m0
			if abs(best[i])>= 0.25 #&& lest[i] >= lmin && lest[i] <=lmax
				mactual=mactual+1
			end
		end
		
		resvecs=A*x-B*x*spdiagm(lest)
		reslist=zeros(m0)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		#(res,ind)=findmax(reslist)
		reslistsort=sort(reslist)
		
		if mactual>0
			res=reslistsort[mactual]		
		else
			res=reslistsort[m0]
		end			

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


function feast_lin_pre(A,B,x0,nc,lmin,lmax,m0,eps,maxit,linsolve)
	
	#linsolve(A,b) ->returns solution to A*x=b

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

		linresnc=zeros(nc)
		linitsnc=zeros(nc)	
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_point[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			G=zk*B-A
			
			Qk=linsolve(G,B*x,x)#\(G,x)

			Q=Q-(gauss_weight[k]/2)*real(r*exp(im*theta)*Qk)
		end

		#=for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end
		=#

		#(Q,R)=qr(Q)

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*B*Q)		
		
		try
			F=eigfact(Hermitian(Aq),Hermitian(Bq))
		catch theerror
			F=eigfact(Bq)
			lest=F[:values]
			println("Bq eigenvalues:\n $(lest)\n")
			error("Error in reduced eigenvalue problem\n Original error:\n $theerror")
		end

		#(val,vec)=eig(Hermitian(Bq))
		#println("   $val")

		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		F=eigfact(Bq)
		best=F[:values]
		mactual=0
		for i in 1:m0
			if abs(best[i])>= 0.25 #&& lest[i] >= lmin && lest[i] <=lmax
				mactual=mactual+1
			end
		end

		#println("    B svs: ",best)
		
		resvecs=A*x-B*x*spdiagm(lest)
		reslist=zeros(m0)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		#(res,ind)=findmax(reslist)
		reslistsort=sort(reslist)
		
		if mactual>0
			res=reslistsort[mactual]		
		else
			res=reslistsort[m0]
		end			

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
