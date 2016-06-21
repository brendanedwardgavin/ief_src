
function myQR(A::Array{Complex128,2})

        (n,m)=size(A)

        Q=zeros(Complex128,n,m)

        Q[:,1]=A[:,1]
	Q[:,1]=Q[:,1]/norm(Q[:,1])

        for i in 2:m
                Q[:,i]=A[:,i]
                for j in 1:i
                        Q[:,i]=Q[:,i]-(Q[:,j]'*Q[:,i])[1]*Q[:,j]
                end
                Q[:,i]=Q[:,i]/norm(Q[:,i])
        end

        R=Q'*A

        return Q,R

end

function powerSolve(A,x0,its)
	x=x0	

	for i in 1:its
		x=A*x
		(q,r)=qr(x)
		x=q	
	end

	return x
end

function improveSolveSub(A1::Array{Complex128,2},B1::Array{Complex128,2},r,m0,eps,Pi)
	(n,m)=size(B1)

	its2=0
	X=zeros(n,m)
	residual=0.0
	for i in 1:n
		its2=i
		
		(X,its,residual)=improveSolvePeps(A1,B1,1,i,eps,Pi)
		if residual<eps
			break
		end
	end
	println("    Lin sys res=$residual, its=$its2")
	return X,its2,residual
end

function improveSolvePeps(A1::Array{Complex128,2},B1::Array{Complex128,2},r,m0,eps,Pi)

	(n,m)=size(B1)

	X=zeros(Complex128,n,m)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)	

	Piinv=inv(Pi)
	A=A1*Piinv
	B=B1
	Bc[:]=B
	residual=eps+1.0
	#for i in 1:r		
	its=0
	while residual>eps && its<r
		its=its+1
		
		Xc=lsSolve(A,Bc,m0)	
		X=X+Piinv*Xc
		#X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
		Bc=B-A1*X

		resnorms=zeros(m)
		for j in 1:m
			resnorms[j]=norm(Bc[:,j])/norm(B[:,j])
		end
		residual=maximum(resnorms)
	end

	#println("    Lin sys res=$residual, its=$its")

	its=its*m0
	return X,its,residual
end


function improveSolveP(A1::Array{Complex128,2},B1::Array{Complex128,2},r,m0)

	(n,m)=size(B1)

	X=zeros(Complex128,n,m)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)	

	#Pi=eye(n)
	P=zeros(Complex128,n,n)
	#mid=50
	#P[1:mid,1:mid]=A1[1:mid,1:mid]
	#P[mid:n,mid:n]=A1[mid:n,mid:n]
	P=eye(n)
	P=diagm(diag(A1))
	Pi=A1*B1*B1'+P*(eye(n)-B1*B1')
	#Pi=P
	#Pi=B1*(A1*B1)'+eye(n)
	#Pi=eye(n)
	#Pi=B1*diagm(lest)*(B1)'+eye(n)

	#Piinv=eye(n)-B1*inv(eye(m)+B1'*A1'*B1)*B1'*A1'
	Piinv=inv(Pi)
	#Piinv=eye(n)
	#Pi=A1*B1*B1'+eye(n)
	#Pi=eye(n)
	#Pi=P
	A=A1*Piinv
	#A=A1*inv(Pi)

	#A=inv(Pi)*A1
	#B=inv(Pi)*B1
	B=B1
	#A=A1
	#Pi=eye(n)
	for i in 1:r
		
		#Bc=B-A*X
		Bc=B-A1*X

		#println("$i: norm bc=$(norm(Bc))")

		Xc=lsSolve(A,Bc,1)	
		X=X+Piinv*Xc
		#X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	return X
end


function improveSolveR(A,B,r,m0,Res)

	V=[B Res]
	lsA=A*V

	(Q,R)=qr(lsA)
	x=V*\(R,Q'*B)

	return x

end

function improveSolve(A,B,r,m0)

	(n,m)=size(B)

	X=zeros(Complex128,n,m)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)
	Bc[:]=B

	for i in 1:r
		Bc=B-A*X

		#println("$i: norm bc=$(norm(Bc))")

		Xc=lsSolve(A,Bc,m0)
		
		rc=Bc-A*Xc

		X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	#=res=zeros(m)
	Bc=B-A*X
	for i in 1:m
		res[i]=norm(Bc[:,i])/norm(B[:,i])	
	end
	println("")
	println("Error = $(maximum(res))")
	println("")
	error("done")=#

	return X
end

function improveSolveIn(A::Array{Complex128,2},B::Array{Complex128,2},x0::Array{Float64,2},r,m0)

	(n,m)=size(B)

	#X=zeros(Complex128,n,m)
	X=x0
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)
	Bc[:]=B

	for i in 1:r
		Bc=B-A*X

		#println("$i: norm bc=$(norm(Bc))")

		Xc=lsSolve(A,Bc,m0)
		
		rc=Bc-A*Xc

		X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	return X
end

function lsSolve(A::Array{Complex128,2},B::Array{Complex128,2},m0)

	
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
		Q[:,m*k+1:m*(k+1)]=A*Q[:,m*(k-1)+1:m*k]
		for i in k*m+1:(k+1)*m
                	Q[:,i]=Q[:,i]/norm(Q[:,i])
			#=for j in 1:i-1
                        	c=(Q[:,i]'*Q[:,j])[1]
	                        Q[:,i]=Q[:,i]-c*Q[:,j]
        	                Q[:,i]=Q[:,i]/norm(Q[:,i])
                	end=#
        	end
	end	
	As=A*Q
	
	(W,R)=qr(As)
	#(W,R)=myQR(As)
	
	Bs=W'*B
	Xs=\(R,Bs)
	X=Q*Xs
	return X
	#X=Q*Xs
end

function improveSolveG(A::Array{Complex128,2},Bm::Array{Float64,2},B::Array{Complex128,2},r,m0)

	(n,m)=size(B)

	invB=inv(Bm)
	#invB=eye(n)

	X=zeros(Complex128,n,m)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)
	Bc[:]=B

	for i in 1:r
		Bc=B-invB*A*X
		#Bc=Bm*B-A*X

		Xc=lsSolveG(A,Bm,Bc,m0,X)
		
		#rc=Bc-A*Xc

		X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	return X
end

function cgB(A,b,x0,it,eps)
	x=x0

	r=b-A*x
	p=r
	rsold=(r'*r)[1]
		
	nits=0
	for k in 1:it
		Ap=A*p
		alpha=rsold/(p'*Ap)[1]
		x=x+alpha*p
		r=r-alpha*Ap
		rsnew=(r'*r)[1]
		
		#println("$k: $(sqrt(rsnew))")
		if sqrt(abs(rsnew))<eps
			break
		end
		p=r+(rsnew/rsold)*p
		rsold=rsnew
		nits=k
	end
	#println(nits)
	return x
end

function improveSolveGIn(A::Array{Complex128,2},Bm::Array{Float64,2},B::Array{Complex128,2},x0::Array{Float64,2},r,m0)

	(n,m)=size(B)

	invB=inv(Bm)
	#invB=eye(n)

	#X=zeros(Complex128,n,m)
	X=convert(Array{Complex128,2},x0)
	Xc=zeros(Complex128,n,m)
	Bc=zeros(Complex128,n,m)
	Bc[:]=B

	for i in 1:r
		ax=A*X
		if i>1
			startx=B
		else
			startx=rand(n,m)
		end
		ibax=ax
		for j in 1:m
		#	ibax[:,j]=cgB(Bm,ax[:,j],startx[:,j],100,1e-10)
		end		

		

		#Bc=B-ibax
		#Bc=B-invB*A*X
		Bc=Bm*B-A*X

		Xc=lsSolveG(A,Bm,Bc,m0,X)
		
		#rc=Bc-A*Xc

		X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	return X
end

function bnrm(v::Array{Complex128,1},B::Array{Float64,2})
        return sqrt((v'*B*v)[1])
end

function myGQR(A::Array{Complex128,2},B::Array{Float64,2})

        (n,m)=size(A)

        Q=zeros(Complex128,n,m)

        Q[:,1]=A[:,1]/bnrm(A[:,1],B)

        for i in 2:m
                Q[:,i]=A[:,i]
                for j in 1:i
                        Q[:,i]=Q[:,i]-(Q[:,i]'*B*Q[:,j])[1]*Q[:,j]/bnrm(Q[:,j],B)
                end
                Q[:,i]=Q[:,i]/bnrm(Q[:,i],B)
        end

        R=Q'*B*A

        return Q,R

end

function lsSolveG(A::Array{Complex128,2},Bm::Array{Float64,2},B::Array{Complex128,2},m0,Xi::Array{Complex128,2})
	
	
	(n,m)=size(B)
	invB=inv(Bm)
	#invB=eye(n)
	X=zeros(Complex128,n,m)	

	Q=zeros(Complex128,n,m*m0)
	Q[:,1:m]=B

        Q[:,1]=Q[:,1]/norm(Q[:,1])
        for i in 2:m	
                #Q[:,i]=Q[:,i]/norm(Q[:,i])
		#=for j in 1:i-1
			c=(Q[:,i]'*Q[:,j])[1]
			Q[:,i]=Q[:,i]-c*Q[:,j]
			Q[:,i]=Q[:,i]/norm(Q[:,i])
                end=#
        end

	for k in 1:m0-1
		Q[:,m*k+1:m*(k+1)]=inv(Bm)*A*Q[:,m*(k-1)+1:m*k]
		for i in k*m+1:(k+1)*m
                	#Q[:,i]=Q[:,i]/norm(Q[:,i])
			#=for j in 1:i-1
                        	c=(Q[:,i]'*Q[:,j])[1]
	                        Q[:,i]=Q[:,i]-c*Q[:,j]
        	                Q[:,i]=Q[:,i]/norm(Q[:,i])
                	end=#
        	end
	end	

	#B-orthogonalize Q with respect to Xi
	#=if sum(abs(Xi))>0
	for i in 1:m*m0
		for j in 1:m
			Q[:,i]=Q[:,i]-(Q[:,i]'*Bm*Xi[:,j])[1]*Xi[:,j]/bnrm(Xi[:,j],Bm)
		end
		Q[:,i]=Q[:,i]/bnrm(Q[:,i],Bm)	
	end
	end=#

	if isnan(sum(Q))
		error("Q is NAN")
	end

	startx=zeros(n)
	for j in 1:m
		Q[:,j]=cgB(Bm,Q[:,j],startx,100,1e-13)
	end
	
	As=A*Q	

	(W,R)=qr(As)
	#(W,R)=myGQR(As,Bm)	
	
	Bs=W'*B
	
	#Bs=W'*Bm*B
	Xs=\(R,Bs)
	X=Q*Xs
	return X
	#X=Q*Xs
end
