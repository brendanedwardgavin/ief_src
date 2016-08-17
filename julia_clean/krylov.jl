function basic_krylov(A,X,k)

	(n,m)=size(x0)
	
	V=zeros(n,m*k)
	V[:,1:m]=X
	(qm,rm)=qr(V[:,1:m])
	V[:,1:m]=qm
	for i in 2:k
		Vold=V[:,1:(i-1)*m]
		Vnew=A*V[:,1+(i-2)*m:(i-1)*m]
		psi=\(Vold'*Vold,Vold'*Vnew)
		Vnew=Vnew-Vold*psi
		(qm,rm)=qr(Vnew)
		V[:,1+(i-1)*m:i*m]=qm
	end

	#(Q,R)=qr(V)
	#V[:]=Q

	Av=V'*A*V
	Bv=V'*V

	(lest,xest)=eig(Av,Bv)

	Xnew=V*xest
	lestout=lest
	#res=A*Xnew-Xnew*diagm(lest)


	return (lestout,Xnew)
end


function folded_krylov(A1,X,k,emin,emax)
	(n,m)=size(x0)

	emid=((emax+emin)/2.0)
	A=(A1-emid*eye(n))^2	

	
	V=zeros(n,m*k)
	V[:,1:m]=X
	(qm,rm)=qr(V[:,1:m])
	V[:,1:m]=qm
	for i in 2:k
		Vold=V[:,1:(i-1)*m]
		Vnew=A*V[:,1+(i-2)*m:(i-1)*m]
		psi=\(Vold'*Vold,Vold'*Vnew)
		Vnew=Vnew-Vold*psi
		(qm,rm)=qr(Vnew)
		V[:,1+(i-1)*m:i*m]=qm
	end

	#(Q,R)=qr(V)
	#V[:]=Q

	Av=V'*A*V
	Bv=V'*V

	(lest,xest)=eig(Av,Bv)

	Xnew=V*xest

	rx=A*Xnew-Xnew*diagm(lest)
	resE=zeros(m*k,2)
	for i in 1:m*k
	        resE[i,1]=norm(rx[:,i])/norm(Xnew[:,i])
        	resE[i,2]=i
	end

	sortedResE=sortrows(resE,by=x->x[1])
	
	xsel=zeros(n,m*k)
	for i in 1:m*k
		xsel[:,i]=Xnew[:,convert(Int64,sortedResE[i,2])]
	end

	Av=xsel'*A1*xsel
	Bv=xsel'*xsel
	(lout,xselsm)=eig(Av,Bv)

	Xout=xsel*xselsm
	return (lout,Xout)
end


function restart_krylov(A,X,k,restarts)

	(n,m)=size(x0)
	Xnew=zeros(n,k*m)
	Xnew[:,1:m]=X
	lest=zeros	

	V=zeros(n,m*k)
	for j in 1:restarts
		V[:,1:m]=Xnew[:,1:m]
		(qm,rm)=qr(V[:,1:m])
		V[:,1:m]=qm
		for i in 2:k
			Vold=V[:,1:(i-1)*m]
			Vnew=A*V[:,1+(i-2)*m:(i-1)*m]
			psi=\(Vold'*Vold,Vold'*Vnew)
			Vnew=Vnew-Vold*psi
			(qm,rm)=qr(Vnew)
			V[:,1+(i-1)*m:i*m]=qm
		end

		#(Q,R)=qr(V)
		#V[:]=Q

		Av=V'*A*V
		Bv=V'*V

		(lest,xest)=eig(Av,Bv)

		Xnew=V*xest
	end

	#res=A*Xnew-Xnew*diagm(lest)

	return (lest,Xnew)
end



function mid_krylov(A,X0,k)

	(n,m)=size(x0)


	ksize=20
	V=zeros(n,m*ksize)
	
	X=zeros(n,m)
	X[:]=X0

	lest=zeros(1:m*ksize)
	lestout=zeros(m)

	for j in 1:k
		println("Outer iteration $j of $k")
		V[:,1:m]=X
		(qm,rm)=qr(V[:,1:m])
		V[:,1:m]=qm
		for i in 2:ksize
			#println("    Inner iteration $i of $ksize")
			Vold=V[:,1:(i-1)*m]
			Vnew=A*V[:,1+(i-2)*m:(i-1)*m]
			psi=\(Vold'*Vold,Vold'*Vnew)
			Vnew=Vnew-Vold*psi
			(qm,rm)=qr(Vnew)
			V[:,1+(i-1)*m:i*m]=qm
		end
		
		#(Q,R)=qr(V)
		#V[:]=Q

		Av=V'*A*V
		Bv=V'*V

		(lest,xest)=eig(Av,Bv)

		sortmat=zeros(1+m*ksize,m*ksize)
		sortmat[1,1:m*ksize]=lest
		sortmat[2:m*ksize+1,1:m*ksize]=xest
		sorted=sortcols(sortmat,by=x->x[1])
	
		lowindex=floor(Int,m*ksize/2-m/2)
		highindex=lowindex+m-1
		midindex=floor(Int,m*ksize/2)
		#X=V*sorted[2:m*ksize+1,lowindex:highindex]
		#lestout=lest[lowindex:highindex]
		X=V*sorted[2:m*ksize+1,lowindex:highindex]
		lestout=vec(sorted[1,lowindex:highindex])
		#println("midindex = $midindex, lowindex = $lowindex, highindex = $highindex")
		#error("blerg")
	end

	#res=A*Xnew-Xnew*diagm(lest)

	return (lestout,X)
end
