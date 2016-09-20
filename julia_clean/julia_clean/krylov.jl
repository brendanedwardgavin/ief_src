function fomfeastBlock(A,x0,lmin,lmax,eps,maxit,ksize)

	(n,m)=size(x0)

	X=zeros(n,m)
	X[:]=x0

	eigsret=zeros(m)

	for i in 1:maxit

		#for each vec:
		
		#form subspace
		Q=zeros(n,(ksize+1)*m)
		Q[:,1:m]=X[:,1:m]
	
		for k in 2:ksize+1
			Q[:,(k-1)*m+1:k*m]=A*Q[:,(k-2)*m+1:(k-1)*m]
		end

		(V,R)=qr(Q)
		
	#solve reduced problem and project

		Av=V'*A*V
		(eigvals,vecs)=eig(Av)
		neigs=0
		for k in 1:(ksize+1)*m
			if eigvals[k]<lmin || eigvals[k]>lmax
				eigvals[k]=0.0
			else
				eigvals[k]=1.0
				neigs=neigs+1
			end
		end

		println("    ---$neigs")
		
		
		P=vecs*diagm(eigvals)*vecs'

		X=V*P*V'*X
		
		
		#solve reduced problem...again, and check error
		(Q,R)=qr(X)
		X[:]=Q	
		Av=X'*A*X
		Bv=X'*X

		(eigsret,vecs)=eig(Av,Bv)

		X=X*vecs

		res=A*X-X*diagm(eigsret)
		err=norm(res)/norm(X)
		println("fomfeast $i      $err")
		println("    $(sort(eigsret))")
	end

	return (eigsret,X)
	
end



function gmresfeast(A,x0,lmin,lmax,eps,maxit,ksize)

	(n,m)=size(x0)

	X=zeros(n,m)
	X[:]=x0

	eigsret=zeros(m)

	emid=(lmin+lmax)/2.0

	for i in 1:maxit

		#for each vec:
		for j in 1:m
		#form subspace
			Q=zeros(n,ksize+1)
			Q[:,1]=X[:,j]
		
			for k in 2:ksize+1
				Q[:,k]=A*Q[:,k-1]
			end

			(V,R)=qr(Q)
			(W,R)=qr(A*Q)
		#solve reduced problem and project
			

			#Av=V'*A*V
			Av1=W'*A*V
			S=W'*V
			Y=\(S',Av1')
			
			Av=Y'
			(eigvals,vecs)=eig(Av)

			neigs=0	
			for k in 1:ksize+1
				#if eigvals[k]<lmin || eigvals[k]>lmax
				if abs(eigvals[k]-emid)>0.5*abs(lmax-lmin)
					eigvals[k]=0.01
				else
					eigvals[k]=1.0
					neigs=neigs+1
				end
			end

			println("    ---$neigs")			

			P=vecs*diagm(eigvals)*inv(vecs)

			T=\(S,W'*X[:,j])			

			X[:,j]=real(V*P*T)	
			
		end

		#solve reduced problem...again, and check error
		(Q,R)=qr(X)
		X[:]=Q	
		Av=X'*A*X
		Bv=X'*X

		(eigsret,vecs)=eig(Av,Bv)

		X=X*vecs

		res=A*X-X*diagm(eigsret)
		err=norm(res)/norm(X)
		println("fomfeast $i      $err")
		println("    $(sort(eigsret))")
	end

	return (eigsret,X)
	
end




function fomfeast(A,x0,lmin,lmax,eps,maxit,ksize)

	(n,m)=size(x0)

	X=zeros(n,m)
	X[:]=x0

	eigsret=zeros(m)

	emid=(lmin+lmax)/2.0

	for i in 1:maxit

		#for each vec:
		for j in 1:m
		#form subspace
			Q=zeros(n,ksize+1)
			Q[:,1]=X[:,j]
		
			for k in 2:ksize+1
				Q[:,k]=A*Q[:,k-1]
			end

			(V,R)=qr(Q)
		#solve reduced problem and project
			
			Av=V'*A*V
			#Av=W'*A*V
			(eigvals,vecs)=eig(Av)
			neigs=0
			for k in 1:ksize+1
				#if eigvals[k]<lmin || eigvals[k]>lmax
				if abs(eigvals[k]-emid)>0.5*abs(lmax-lmin)
					#eigvals[k]=0.001
				else
					#eigvals[k]=1.0
					neigs=neigs+1
				end

				eigvals[k]=1/abs(eigvals[k]-emid)
			end

			println("    ---$neigs")
			
			
			P=vecs*diagm(eigvals)*vecs'
			#P=vecs*diagm(eigvals)*

			X[:,j]=V*P*V'*X[:,j]	
			
			
		end

		#solve reduced problem...again, and check error
		(Q,R)=qr(X)
		X[:]=Q	
		Av=X'*A*X
		Bv=X'*X

		(eigsret,vecs)=eig(Av,Bv)

		X=X*vecs

		res=A*X-X*diagm(eigsret)
		err=norm(res)/norm(X)
		println("fomfeast $i      $err")
		println("    $(sort(eigsret))")
	end

	return (eigsret,X)
	
end


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
