function psolve(A,x0,eps,maxit)

	(n,m)=size(x0)

	residuals=zeros(maxit)

	x=x0
	(q,r)=qr(x)
	x=q
	
	init=5

	lest=zeros(m)

	for i in 1:maxit

		for j in 1:m
			x[:,j]=x[:,j]/norm(x[:,j])
		end	
		As=x'*A*x
		Bs=x'*x
		F=eigfact(As,Bs)
                lest=F[:values]
                xest=F[:vectors]
                x=x*xest
		
		xb=x*diagm(lest)
		#x=zeros(n,m)
		for j in 1:init
			#W=zeros(n,2*m)
			#W[:,1:m]=xb-A*x
			#W[:,m+1:2*m]=A*(xb-A*x)
			W=xb-A*x
			Aw=A*W
			(Q,R)=qr(Aw)
			
			ws=Q'*W[:,1:m]
			xw=\(R,ws)
			x=x+W*xw
		end

		resvecs=A*x-x*diagm(lest)
                reslist=zeros(m,1)
                for j in 1:m
                        reslist[j]=norm(resvecs[:,j])/norm(x[:,j])
                end
                (res,ind)=findmin(reslist)
		residuals[i]=res
		#println("$i   $res   $(lest[1])")	
		if res<eps
			break
		end
	end

	return (lest,x,residuals)
end


function psolveInv(A,x0,eps,maxit)

	(n,m)=size(x0)

	residuals=zeros(maxit)

	x=x0
	(q,r)=qr(x)
	x=q
	
	init=5

	lest=zeros(m)

	for i in 1:maxit

		for j in 1:m
			x[:,j]=x[:,j]/norm(x[:,j])
		end	
		As=x'*A*x
		Bs=x'*x
		F=eigfact(As,Bs)
                lest=F[:values]
                xest=F[:vectors]
                x=x*xest
		
		xb=x*inv(diagm(lest))
		#x=zeros(n,m)
	
		#calculate dx:
		dx=A*xb-x

		x=x+dx	

		resvecs=A*x-x*diagm(lest)
                reslist=zeros(m,1)
                for j in 1:m
                        reslist[j]=norm(resvecs[:,j])/norm(x[:,j])
                end
                (res,ind)=findmin(reslist)
		residuals[i]=res
		#println("$i   $res   $(lest[1])")	
		if res<eps
			break
		end
	end

	return (lest,x,residuals)
end
