function blockGivens(vec)

	(m,n)=size(vec)
	vectype=typeof(vec[1])

	M=zeros(2*n,2*n)

	Y=vec[n+1:2*n,1:n]
	X=vec[1:n,1:n]

	B=\(X',Y')'

	(l,x)=eig(B'*B)
	for i in 1:n
		l[i]=sqrt(1/(1+l[i]))
	end
	A=diagm(l)*x'

	(l,x)=eig(B*B')
	for i in 1:n
		l[i]=sqrt(1/(1+l[i]))
	end
	C=diagm(l)*x'

	M[1:n,1:n]=A
	M[1:n,n+1:2*n]=A*B'
	M[n+1:2*n,1:n]=-1*C*B
	M[n+1:2*n,n+1:2*n]=C

	return M
end

function gmres2(A,B,k,r,eps)

	(n,m)=size(B)
	mattype=typeof(A[1,1])
	

	R=zeros(mattype,n,m)
	X=zeros(mattype,n,m)
	X0=zeros(mattype,n,m)
	
	V=zeros(mattype,n,m*(k+1))
	H=zeros(mattype,m*(k+1),m*k)
	Bsm=zeros(mattype,m*(k+1),m)

	err=0.0
	for i in 1:r
		X0[:]=X
		R=B-A*X		

		V[:,1:m]=R

		println("V=",V[1,1])
		for j in 1:k	
			#blockArnoldiPlus(A,V,H,j,k,m)
			blockArnoldiIt(A,V,H,Bsm,j,k,m)
		
			println("H=",sum(H))

			println("Bsm=",sum(Bsm))

			println("V=",sum(V))
	
			ym=\(H[1:(j+1)*m,1:j*m],Bsm[1:(j+1)*m,1:m])

			#(Q,R)=qr(H[1:(j+1)*m,1:j*m])
			#ym=\(R,Q'*Bsm[1:(j+1)*m,1:m])

			println("Ym=",ym[1,1])

			X=X0+V[1:n,1:j*m]*ym
			#Res=B-A*X
			Res=B-A*X0-V[1:n,1:(j+1)*m]*H[1:(j+1)*m,1:j*m]*ym
			
			println(i," ",j)
			for l in 1:m
				println("     ", norm(Res[:,l])/norm(B[:,l]))
			end
			err=norm(Res)/norm(B)
			#println("$i, $j: $err")
			if err<eps
				break
			end
		end
		
		if err<eps
			break
		end
	end
	
	return X
end



function gmres3(A,B,k,r,eps)
#GMRES by doing full krylov subspace at each restart
	(n,m)=size(B)
	mattype=typeof(A[1,1])
	

	R=zeros(mattype,n,m)
	X=zeros(mattype,n,m)

	for i in 1:r
		R=B-A*X

		(V,H,Bsm)=blockArnoldi(A,R,k)
		ym=\(H,Bsm)
		X=X+V[1:n,1:m*k]*ym

		Res=B-A*X
		err=norm(Res)/norm(B)
		println("$i: ",err)
		if err<eps
			break
		end
	end
	
	return X
end




function blockArnoldiStart(A,V0,k)
#set everything up and do the first iteration of arnoldi
	(n,m)=size(V0)
	Vtype=typeof(V0[1,1])
	V=zeros(Vtype,n,m*(k+1))
	H=zeros(Vtype,m*(k+1),m*k)

	Bsm=zeros(Vtype,m*(k+1),m)
	(Q,R)=qr(V0)
	Bsm[1:m,1:m]=R

	V[:,1:m]=Q
	
	i=1
	
	Vnew=A*V[:,(i-1)*m+1:i*m]
	for j in 1:i
		Hnew=V[:,(j-1)*m+1:j*m]'*Vnew
		H[(j-1)*m+1:j*m,(i-1)*m+1:i*m]=Hnew

		Vnew=Vnew-V[:,(j-1)*m+1:j*m]*Hnew
	end
	(Q,R)=qr(Vnew)
	V[:,i*m+1:(i+1)*m]=Q
	
	H[i*m+1:(i+1)*m,(i-1)*m+1:i*m]=R


	return (V,H,Bsm)

end



function blockArnoldiPlus(A,V0,H0,k0,k,m)
#assume V0, H0 are already allocated correctly
#do the next iteration of arnoldi
	n=size(V0,1)
	#m=convert(Int64,m/k)

	Vtype=typeof(V0[1,1])

	i=k0

	Vnew=A*V0[:,(i-1)*m+1:i*m]
	for j in 1:i
		Hnew=V0[:,(j-1)*m+1:j*m]'*Vnew
		H0[(j-1)*m+1:j*m,(i-1)*m+1:i*m]=Hnew

		Vnew=Vnew-V0[:,(j-1)*m+1:j*m]*Hnew
	end
	(Q,R)=qr(Vnew)
	V0[:,i*m+1:(i+1)*m]=Q
	H0[i*m+1:(i+1)*m,(i-1)*m+1:i*m]=R

end


function blockArnoldiIt(A,V,H,Bsm,k0,k,m)
#assume V0, H0 are already allocated correctly
#do the next iteration of arnoldi
	n=size(V,1)
	#m=convert(Int64,m/k)

	Vtype=typeof(V[1,1])

	if k0==1	
		(Q,R)=qr(V[:,1:m])

		V[:,1:m]=Q
		Bsm[1:m,1:m]=R[1:m,1:m]
	end



	i=k0

	Vnew=A*V[:,(i-1)*m+1:i*m]
	for j in 1:i
		Hnew=V[:,(j-1)*m+1:j*m]'*Vnew
		H[(j-1)*m+1:j*m,(i-1)*m+1:i*m]=Hnew

		Vnew=Vnew-V[:,(j-1)*m+1:j*m]*Hnew
	end
	(Q,R)=qr(Vnew)
	V[:,i*m+1:(i+1)*m]=Q
	H[i*m+1:(i+1)*m,(i-1)*m+1:i*m]=R

end



function blockArnoldi(A,V0,k)

	(n,m)=size(V0)
	Vtype=typeof(V0[1,1])
	V=zeros(Vtype,n,m*(k+1))
	H=zeros(Vtype,m*(k+1),m*k)

	(Q,R)=qr(V0)

	Bsm=zeros((k+1)*m,m)
	Bsm[1:m,1:m]=R

	V[:,1:m]=Q
	for i in 1:k
		Vnew=A*V[:,(i-1)*m+1:i*m]
		for j in 1:i
			Hnew=V[:,(j-1)*m+1:j*m]'*Vnew
			H[(j-1)*m+1:j*m,(i-1)*m+1:i*m]=Hnew

			Vnew=Vnew-V[:,(j-1)*m+1:j*m]*Hnew
		end
		(Q,R)=qr(Vnew)
		V[:,i*m+1:(i+1)*m]=Q
		H[i*m+1:(i+1)*m,(i-1)*m+1:i*m]=R
	end


	return (V,H,Bsm)
end




function blockArnoldiGivens(A,V0,k)

	(n,m)=size(V0)
	Vtype=typeof(V0[1,1])
	V=zeros(Vtype,n,m*(k+1))
	H=zeros(Vtype,m*(k+1),m*k)

	(Q,R)=qr(V0)

	Bsm=zeros((k+1)*m,m)
	Bsm[1:m,1:m]=R

	V[:,1:m]=Q
	for i in 1:k
		Vnew=A*V[:,(i-1)*m+1:i*m]
		for j in 1:i
			Hnew=V[:,(j-1)*m+1:j*m]'*Vnew
			H[(j-1)*m+1:j*m,(i-1)*m+1:i*m]=Hnew

			Vnew=Vnew-V[:,(j-1)*m+1:j*m]*Hnew
		end
		(Q,R)=qr(Vnew)
		V[:,i*m+1:(i+1)*m]=Q
		H[i*m+1:(i+1)*m,(i-1)*m+1:i*m]=R
	end


	return (V,H,Bsm)
end





function affine(A,B,maxit)

	(n,m)=size(B)

	X=zeros(n,m)
	X[:]=B

	block=3

	for i in 1:maxit
		#k=rand(1:n)
		#beta=(B[k,1]-(A[k,:]*X)[1])/(A[k,:]*A[k,:]')[1]
		#X=X+A[k,:]'*beta		
	
		ks=randperm(n)
		rowvecs=zeros(n,block)
		Bs=zeros(block,m)
		for j in 1:block
			rowvecs[:,j]=A[ks[j],:]'
			Bs[j,1]=B[ks[j],1]
		end
		beta=\(rowvecs'*rowvecs, Bs-rowvecs'*X)
		X=X+rowvecs*beta
		
		err=norm(A*X-B)/norm(B)

		println("Affine it $i    $err")
	end

	return X
end



function fom(A,B,maxit,eps)

	(n,m)=size(B)

	V=zeros(Complex128,n,maxit+1)

	X=zeros(Complex128,n,m)

	for i in 1:m
		V[:,1]=B[:,1]/norm(B[:,1])
		xi=zeros(Complex128,n,1)
		println("----$i")
		for j in 2:maxit+1
			Vj=A*V[:,j-1]	
			V[:,j]=Vj

			(Q,R)=qr(V[:,1:j])
			#Q=V[:,1:j]
			#xi=Vj*\(Vj'*A*Vj,Vj'*B[:,i])
			#xi=V[:,1:j]*\(V[:,1:j]'*A'*A*V[:,1:j],V[:,1:j]'*A'*B[:,i])
			xi=Q*\(Q'*A*Q,Q'*B[:,i])

			err=norm(B[:,i]-A*xi)/norm(B[:,i])
			if err<eps
				 break
			end
			println(err)
		end

		X[:,i]=xi
	end
	
	return X
end


function blockCGbasic(A1,B1,its)

	A=A1#A1'*A1
	B=B1#A1'*B1

	(n,m)=size(B)
	
	X=zeros(n,m)
	R=A*X-B
	P=-1.0*R
	
	for i in 1:its
		lambda=inv(P'*A*P)*R'*R
		X=X+P*lambda
		Rnew=R+A*P*lambda
		psi=inv(R'*R)*Rnew'*Rnew
		P=-1.0*Rnew+P*psi
		R[:]=Rnew

		#println("$i  res=",norm(B1-A1*X)/norm(B1))
	end	

	return X
end



function blockCGbasicnorm(A,B1,its)

	B=A'*B1

	(n,m)=size(B)

	X=zeros(n,m)
	R=A'*(A*X)-B
	P=-1.0*R

	#println("X=",X,"\nR=",R,"\nP=",P,"\n\n")
	
	for i in 1:its
		lambda=\((A*P)'*(A*P),R'*R)
		X=X+P*lambda
		Rnew=R+A'*(A*P)*lambda
		psi=\(R'*R,Rnew'*Rnew)
		P=-1.0*Rnew+P*psi
		R[:]=Rnew
		
		#println("Lambda=",lambda,"\nXnew=",X,"\nRnew=",R,"\npsi=",psi,"\nPnew=",P)
	
		#println("\n\n$i  res=",norm(B1-A*X)/norm(B1),"\n\n")
	end	

	return X
end

function blockBCG(A,B,M,its)

	(n,m)=size(B)
	
	gamma1=eye(m,m)
	gamma2=eye(m,m)

	X=zeros(Complex128,n,m)

	R1=B-A*X
	R2=B-A'*X
	P1=M*R1*gamma1
	P2=M'*R2*gamma2

	for i in 1:its
		
		alpha1=inv(P2'*A*P1)*gamma2'*R2'*M*R1

		alpha2=inv(P1'*A'*P2)*gamma1'*R1'*M'*R2

		X=X+P1*alpha1

		R1new=R1-A*P1*alpha1
		R2new=R2-A'*P2*alpha2

		beta1=inv(gamma1)*inv(R2'*M*R1)*R2new'*M*R1new
		beta2=inv(gamma2)*inv(R1'*M'*R2)*R1new'*M'*R2new

		P1=(M*R1new+P1*beta1)*gamma1
		P2=(M'*R2new+P2*beta2)*gamma2
	
		R1[:]=R1new
		R2[:]=R2new

		#println("$i   $(norm(R1)/norm(B))")	
	end

	R=B-A*X
	Xc=gmres_lsSolve(A,B,eye(n),1)
	X=X+Xc		

	return X
end


function blockCGnorm(A,B,M,its)
	
	(n,m)=size(B)
	X=zeros(Complex128,n,m)

	gamma=eye(m,m)

	R=A'*B-A'*A*X
	P=M*R*gamma

	for i in 1:its
		alpha=\((P'*A')*(A*P),gamma'*R'*M*R)
		
		X=X+P*alpha

		newR=R-A'*(A*P)*alpha

		beta=inv(gamma)*inv(R'*M*R)*newR'*M*newR

		P=(M*newR+P*beta)*gamma

		R[:]=newR

		#println("$i   $(norm(R)/norm(A'*B))")
	end

	return X
end


function blockCG(A,B,M,its)
	
	(n,m)=size(B)
	X=zeros(Complex128,n,m)

	gamma=eye(m,m)

	R=B-A*X
	P=M*R*gamma

	for i in 1:its
		alpha=\((P')*(A*P),gamma'*R'*M*R)
		
		X=X+P*alpha

		newR=R-(A*P)*alpha

		beta=inv(gamma)*inv(R'*M*R)*newR'*M*newR

		P=(M*newR+P*beta)*gamma

		R[:]=newR

		#println("$i   $(norm(R)/norm(A'*B))")
	end

	return X
end



function regularCGnormal(A,B,its)
	(n,m)=size(B)
	X=zeros(Complex128,n,m)

	for j in 1:m
		r=A'*B[:,j]-A'*A*X[:,j]
		d=r
		deltanew=(r'*r)[1]
		delta=deltanew

		for i in 1:its
			q=A'*A*d
			alpha=deltanew/(d'*q)[1]
			X[:,j]=X[:,j]+alpha*d

			r=A'*B[:,j]-A'*A*X[:,j]

			delta=deltanew
			deltanew=(r'*r)[1]
			beta=deltanew/delta
			d=r+beta*d
		end

	end

	return X
end


function regularCG(A,B,its)
	(n,m)=size(B)
	X=zeros(Complex128,n,m)

	for j in 1:m
		r=B[:,j]-A*X[:,j]
		d=r
		deltanew=(r'*r)[1]
		delta=deltanew

		for i in 1:its
			q=A*d
			alpha=deltanew/(d'*q)[1]
			X[:,j]=X[:,j]+alpha*d

			r=B[:,j]-A*X[:,j]

			delta=deltanew
			deltanew=(r'*r)[1]
			beta=deltanew/delta
			d=r+beta*d
		end

	end

	return X
end


function gmres(A,B,P,r,m0)
	#restarted GMRES with preconditioner P, matrix A, right hand sides B, krylov size m0, restarts r

	(n,m)=size(B)

	X=zeros(Complex128,n,m)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)
	Bc[:]=B

	Pinv=inv(P)
	A1=A*Pinv

	for i in 1:r
		Bc=B-A*X

		#println("$i: norm bc=$(norm(Bc))")

		Xc=gmres_lsSolve(A,Bc,Pinv,m0)
		
		#rc=Bc-A*Xc

		X=X+Pinv*Xc

	end

	return X
end


function gmres_lsSolve(A,B,Pinv,m0)
	#inner GMRES iteration with inv(preconditioner) Pinv, matrix A, rhs B
	
	(n,m)=size(B)
	X=zeros(Complex128,n,m)	

	Q=zeros(Complex128,n,m*m0)
	Q[:,1:m]=B

        Q[:,1]=Q[:,1]/norm(Q[:,1])
        for i in 2:m	
                Q[:,i]=Q[:,i]/norm(Q[:,i])
		#=for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,i]=Q[:,i]-c*Q[:,j]
			Q[:,i]=Q[:,i]/norm(Q[:,i])
                end=#
        end

	for k in 1:m0-1
		Q[:,m*k+1:m*(k+1)]=A*Pinv*Q[:,m*(k-1)+1:m*k]
		for i in k*m+1:(k+1)*m
                	Q[:,i]=Q[:,i]/norm(Q[:,i])
			#=for j in 1:i-1
                        	c=(Q[:,i]'*Q[:,j])[1]
	                        Q[:,i]=Q[:,i]-c*Q[:,j]
        	                Q[:,i]=Q[:,i]/norm(Q[:,i])
                	end=#
        	end
	end	
	As=A*Pinv*Q
	
	#println("normalE!")
	(W,R)=qr(As)
	Bs=W'*B

	#R=As'*As
	#Bs=As'*B
	
	
	#Xs=regularCG(R,Bs,1)
	#Pcg=eye(m*m0)
	#@time Xs=blockCG(R,Bs,Pcg,1)
	Xs=\(R,Bs)
	
	X=Q*Xs
	return X
	#X=Q*Xs
end



function stationaryIt(A,B,band,its)
	(n,m)=size(B)
	xmax=rand(n)
	xmax=xmax/norm(xmax)
	lmax=0.0
	for i in 1:100
		xmax=A*xmax
		lmax=norm(xmax)
		xmax=xmax/lmax
	end

	Anorm=zeros(n,n)
	Anorm[:]=A#/(1.05*lmax)

	C=zeros(n,n)
	C[:]=Anorm
	N=zeros(n,n)
	for i in band+1:n
		for j in 1:i-band
			C[i,j]=0.0
			C[j,i]=0.0
			N[i,j]=A[i,j]
			N[j,i]=A[j,i]
		end
	end

	#println("A:\n$Anorm")
	#println("C:\n$C")
	#println("N:\n$N")
	#println("lmax=$lmax")
	#error("stop")

	X=zeros(n,m)
	rold=1.0
	for i in 1:its
		X=\(C,B-N*X)

		R=B-A*X
		rnew=(norm(R)/norm(B))
		println(i,"  ",rnew,"    ",rnew/rold)
		rold=rnew
	end

	return X
end
