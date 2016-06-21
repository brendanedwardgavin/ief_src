
function myQR(A)

	(n,m)=size(A)

	Q=zeros(n,m)
	R=zeros(m,m)

	Q[:,1]=A[:,1]/norm(A[:,1])
	
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

function bnrm(v,B)
	return sqrt((v'*B*v)[1])
end

function myGQR(A,B)

	(n,m)=size(A)

	Q=zeros(n,m)
	R=zeros(m,m)

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

srand(1)
n=10
m=3
A=rand(n,m)

bx=rand(n,n)
(bx,r)=qr(bx)
lambda=diagm(rand(n))
lambda=sqrt(lambda)
F=bx*lambda
B=F*F'

(Q,R)=myGQR(A,B)

println("Error:",norm(A-Q*R)/norm(A))

