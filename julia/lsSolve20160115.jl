
function improveSolve(A::Array{Complex128,2},B::Array{Complex128,2},r,m0)

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

	#=for k in 1:m
		Q[:,k]=Q[:,k]/norm(Q[:,k])
	end

	for k in 1:m0-1
		Q[:,m*k+1:m*(k+1)]=A*Q[:,m*(k-1)+1:m*k]
		
		for i in 1:m
			for j in 1:k-1
				cur=k*m+i
				oth=(j-1)*m+i
				
				c=(Q[:,cur]'*Q[:,oth])[1]
	                        Q[:,cur]=Q[:,cur]-c*Q[:,oth]
        	                Q[:,cur]=Q[:,cur]/norm(Q[:,cur])
			end
		end		
	end=#

	#As=Q'*A*Q
	#Bs=Q'*B
	#Xs=\(As,Bs)

	As=A*Q
	#(U,S,V)=svds(As,nsv=m0)
	#sinv=diagm(1.0./S)	
	#Xs=V*sinv*U'*B

	(W,R)=qr(As)
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
		#Bc=B-invB*A*X
		#Bc=Bm*B-A*X
		AX=A*X
		invBAX=zeros(n,m)
		if i>1
			#invBAX=inv(Bm)*AX
			invBAX=improveSolve(convert(Array{Complex128,2},Bm),AX,1,10)
		end
		#Bc=\(Bm,Bc)

		#println("$i: norm bc=$(norm(Bc))")
		Bc=B-invBAX#improveSolve(convert(Array{Complex128,2},Bm),Bc,2,3)

		Xc=lsSolveG(A,Bm,Bc,m0)
		
		#rc=Bc-A*Xc

		X=X+Xc
		#println("$i: ratio=$(norm(rc)/norm(Bc))")
		#println("   Error=",norm(B-A*X)/norm(B))
	end

	return X
end


function lsSolveG(A::Array{Complex128,2},Bm::Array{Float64,2},B::Array{Complex128,2},m0)
	
	
	(n,m)=size(B)
	invB=inv(Bm)
	#invB=eye(n)
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

	#As=invB*A*Q
	As=A*Q
	#(U,S,V)=svds(As,nsv=m0)
	#sinv=diagm(1.0./S)	
	#Xs=V*sinv*U'*B

	(W,R)=qr(As)
	Bs=W'*Bm*B
	Xs=\(R,Bs)
	X=Q*Xs
	return X
	#X=Q*Xs
end
