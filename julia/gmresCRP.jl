function rotmatC(a::Complex128,b::Complex128)
#a=h_ii, b=h_i+1i 
	c=0.0+im*0.0
	s=0.0+0.0*im

	#=if b==0.0
		c=1.0
		s=0.0
	else
		#=if abs(b)>abs(a)
			temp=a/b
			s=1.0/sqrt(1+abs(temp)^2)
			c=temp*s
		else
			temp=b/a
			c=1.0/sqrt(1.0+abs(temp)^2)
			s=temp*c
		end=#
		s=b/sqrt(abs(a)^2+b^2)
		c=a/sqrt(abs(a)^2+b^2)	
	end=#
	s=b'/sqrt(abs(a)^2+abs(b)^2)
	c=a'/sqrt(abs(a)^2+abs(b)^2)
	return c,s

end


function gmresCRP(A,x,b,Minv,restrt,max_it,tol::Float64)

iter=0
flag=0

bnrm2=norm(b)

#r=\(M,b-A*x)
r=b-A*x
error=norm(r)/bnrm2

n=size(A,1)

m=restrt

V=zeros(Complex128,n,m+1)
H=zeros(Complex128,m+1,m)

cs=zeros(Complex128,m,1)
sn=zeros(Complex128,m,1)
e1=zeros(Complex128,n,1)
e1[1]=1.0

for iter in 1:max_it

	#r=\(M,b-A*x)
	r=b-A*x
	V[:,1]=r/norm(r)
	s=norm(r)*e1

	for i in 1:m

		#w=\(M,A*V[:,i])
		w=(A*Minv*V[:,i])

		for k in 1:i
			H[k,i]=(w'*V[:,k])[1]
			w=w-H[k,i]*V[:,k]
		end
		H[i+1,i]=norm(w)
		V[:,i+1]=w/H[i+1,i]
		#givens rotations:
		for k in 1:i-1
			temp=cs[k]*H[k,i]+sn[k]*H[k+1,i]
			H[k+1,i]=-sn[k]*H[k,i]+cs[k]*H[k+1,i]
			H[k,i]=temp
		end
		(cs[i],sn[i])=rotmatC(H[i,i],H[i+1,i])
		temp=cs[i]'*s[i]
		s[i+1]=-sn[i]'*s[i]
		s[i]=temp
		H[i,i]=cs[i]*H[i,i]+sn[i]*H[i+1,i]
		H[i+1,i]=0.0

		error=abs(s[i+1])/bnrm2

		if error<=tol
			#println(size(H[1:i,1:i]),size(s[1:m]))
			y=\(H[1:i,1:i],s[1:i])
			x=x+V[:,1:i]*y
			break
		end

	end

	if abs(error)<=tol
		break
	end

	y=\(H[1:m,1:m],s[1:m])
	x=x+Minv*V[:,1:m]*y
	#r=\(M,b-A*x)
	r=(b-A*x)
	s[m+1]=norm(r)
	error=s[m+1]/bnrm2

	if abs(error)<=tol
		break
	end
end

#println("error=$(abs(error))")

return x

end

function solveSys(A,B,X0,m0,maxit)
	(n,m)=size(X0)
	#n=size(A,1)
	#m=1
	#Q=zeros(n,m)
	#Q[:]=X0

	X=zeros(n,m)
	X[:]=X0
	Bc=B-A*X

	#Q=[X0 A*X0 A*A*X0 A*A*A*X0]

	Q=zeros(Complex128,n,m0*m)

	Q[:,1:m]=B

	for k in 1:maxit	

	#Q=A*X	
	for i in 2:m0
		Q[:,m*(i-1)+1:m*i]=A*Q[:,m*(i-2)+1:m*(i-1)]
	end

	Q[:,1]=Q[:,1]/norm(Q[:,1])
	for i in 2:m*m0
		Q[:,i]=Q[:,i]/norm(Q[:,i])
		for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,j]=Q[:,j]-c*Q[:,i]
			Q[:,j]=Q[:,j]/norm(Q[:,j])
		end
	end
	
	

	#println(size(Q),size(A))

	As=Q'*A*Q
	Bs=Q'*Bc

	#println("B approx error=",norm(Bc-Q*\(Q'*Q,Bs))/norm(Bc))
	#(U,S,V)=svds(Q,nsv=m0)
	#println("Svds=",S)
	Xs=\(As,Bs)
	
	

	#sol=Q*Xs

	#Xs=\(A*Q,Bc)
	sol=Q*Xs

	#sol=gmresRP(A,X0,Bc,eye(n),10,1,1e-6)

	X=X+sol

	Bc=Bc-A*sol

	Q[:,1:m]=Bc
	
	end
	return X
end

function solveSys3(A,B,X0,m0,maxit)
	(n,m)=size(X0)
	#n=size(A,1)
	#m=1
	#Q=zeros(n,m)
	#Q[:]=X0

	X=zeros(n,m)
	X[:]=X0
	Bc=B-A*X

	#Q=[X0 A*X0 A*A*X0 A*A*A*X0]

	Q=zeros(Complex128,n,m0*m)

	Q[:,1:m]=B
	Q[:,1]=Q[:,1]/norm(Q[:,1])
	for i in 2:m
		Q[:,i]=Q[:,i]/norm(Q[:,i])
		for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,j]=Q[:,j]-c*Q[:,i]
			Q[:,j]=Q[:,j]/norm(Q[:,j])
		end
	end


	for k in 1:maxit	

	
	for i in 1:m0-1
		Q[:,m*i+1:m*(i+1)]=A*Q[:,m*(i-1)+1:m*i]

		for i in i*m+1:(i+1)*m
			Q[:,i]=Q[:,i]/norm(Q[:,i])
			for j in 1:i-1
				c=(Q[:,i]'*Q[:,j])[1]
				Q[:,j]=Q[:,j]-c*Q[:,i]
				Q[:,j]=Q[:,j]/norm(Q[:,j])
			end
		end
	end
	#println(size(Q),size(A))

	As=Q'*A*Q
	Bs=Q'*Bc

	#println("B approx error=",norm(Bc-Q*\(Q'*Q,Bs))/norm(Bc))
	#(U,S,V)=svds(Q,nsv=m0)
	#println("Svds=",S)
	Xs=\(As,Bs)
	
	#sol=Q*Xs

	#Xs=\(A*Q,Bc)
	sol=Q*Xs

	#sol=gmresRP(A,X0,Bc,eye(n),10,1,1e-6)

	X=X+sol

	Bc=Bc-A*sol

	Q[:,1:m]=Bc
	
	end
	return X
end

#=
function solveSys2(A,B,X0,m0,maxit)
	(n,m)=size(X0)
	#n=size(A,1)
	#m=1
	#Q=zeros(n,m)
	#Q[:]=X0

	X=zeros(n,m)
	X[:]=X0
	Bc=B-A*X

	#Q=[X0 A*X0 A*A*X0 A*A*A*X0]
	for k in 1:maxit
	Q=zeros(Complex128,n,m0*m)

	Q[:,1:m]=X
	Q[:,1]=Q[:,1]/norm(Q[:,1])
	for i in 2:m
		Q[:,i]=Q[:,i]/norm(Q[:,i])
		for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,j]=Q[:,j]-c*Q[:,i]
			Q[:,j]=Q[:,j]/norm(Q[:,j])
		end
	end

	
	
	for i in 2:m0-1
	Q
	As=Q[:,1:(i-1)*m]'*A*Q
	Bs=Q'*Bc

	println("B approx error=",norm(Bc-Q*\(Q'*Q,Bs))/norm(Bc))
	(U,S,V)=svds(Q,nsv=m0)
	println("Svds=",S)
	Xs=\(As,Bs)
	
	sol=Q*Xs

	X=X+sol

	#Q=A*X
	
	for i in 2:m0
		Q[:,m*(i-1)+1:m*i]=A*Q[:,m*(i-2)+1:m*(i-1)]
	end

	Q[:,1]=Q[:,1]/norm(Q[:,1])
	for i in 2:m
		Q[:,i]=Q[:,i]/norm(Q[:,i])
		for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,j]=Q[:,j]-c*Q[:,i]
			Q[:,j]=Q[:,j]/norm(Q[:,j])
		end
	end

	#println(size(Q),size(A))

	
	Bc=Bc-A*sol

	Q[:,1:m]=Bc
	
	end
	return X
end=#
