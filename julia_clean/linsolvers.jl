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

	println("X=",X,"\nR=",R,"\nP=",P,"\n\n")
	
	for i in 1:its
		lambda=\((A*P)'*(A*P),R'*R)
		X=X+P*lambda
		Rnew=R+A'*(A*P)*lambda
		psi=\(R'*R,Rnew'*Rnew)
		P=-1.0*Rnew+P*psi
		R[:]=Rnew
		
		println("Lambda=",lambda,"\nXnew=",X,"\nRnew=",R,"\npsi=",psi,"\nPnew=",P)
	
		println("\n\n$i  res=",norm(B1-A*X)/norm(B1),"\n\n")
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

