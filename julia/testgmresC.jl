include("gmres.jl")
include("gmresC.jl")

srand(1)

n=100

#A=rand(Complex128,n,n)
A=rand(n,n)
A=A*A'

ls=diagm(1:n)+im*diagm(1:n)
ls=convert(Array{Complex128,2},ls)
#ls=convert(Array{Float64,2},ls)
(Q,R)=qr(A)
#A=Q*ls*Q'
Ac=Q*ls*Q'

#=x=rand(Complex128,n)
b=rand(Complex128,n)
x[:]=b
M=eye(Complex128,n)=#
x=rand(Float64,2*n)
x=x[1:n]
b=rand(Float64,2*n)
b=b[1:n]
#b[n+1:2*n]=0.0
x[:]=b
M=eye(n)
restrt=99
max_it=1
tol=1e-16

A=[real(Ac) -imag(Ac);imag(Ac) real(Ac)]

xs=gmres(Ac,x,b,M,restrt,max_it,tol)

println(norm(Ac*xs-b)/norm(b))

