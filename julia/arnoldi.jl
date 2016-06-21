srand(1)

n=10

k=3
m=2

A=rand(n,n)
A=A+A'

b=rand(n,m)

Q=zeros(n,(k+1)*m)
H=zeros((k+1)*m,k*m)

for i in 1:m
	Q[:,i]=b[:,i]/norm(b[:,i])
end

#=for i in 1:k
	Q[:,i+1]=A*Q[:,i]
	for j in 1:i
		H[j,i]=(Q[:,j]'*Q[:,i+1])[1]
		Q[:,i+1]=Q[:,i+1]-H[j,i]*Q[:,j]
	end
	H[i+1,i]=norm(Q[:,i])
	Q[:,i+1]=Q[:,i+1]/H[i+1,i]
end=#

for i in 1:k
	Q[:,i*m+1:(i+1)*m]=A*Q[:,(i-1)*m+1:i*m]
	for j in 1:i
		H[j,i]=(Q[:,j]'*Q[:,i+1])[1]
		Q[:,i+1]=Q[:,i+1]-H[j,i]*Q[:,j]
	end
	H[i+1,i]=norm(Q[:,i])
	Q[:,i+1]=Q[:,i+1]/H[i+1,i]
end
