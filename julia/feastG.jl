

module juliaFEAST

export feastG

include("gmresRP.jl")
include("lsSolve.jl")

gauss_weight8=[0.362683783378361;0.362683783378361;0.313706645877887;0.313706645877887;0.222381034453374;0.222381034453374;0.101228536290376;0.101228536290376]
gauss_ab8=[0.183434642495649;-0.183434642495649;0.525532409916328;-0.525532409916328;0.796666477413626;-0.796666477413626;0.960289856497536;-0.960289856497536]



function feastG(A::SparseMatrixCSC{Float64,Int64},B::SparseMatrixCSC{Float64,Int64},x0::Array{Float64,2},lmin::Float64,lmax::Float64,m0::Int64,eps::Float64,maxit::Int64)
	nc::Int64=8
	n::Int64=size(x0,1)
	res::Float64=1.0

	r::Float64=(lmax-lmin)/2.0
	x=x0;

	lest=zeros(Float64,m0,1)
	lest=diag(x0'*(A*x0))

	it=0


	maxlinit=3000

	residuals=zeros(Float64,maxit)

	while res>eps && it<maxit
		it=it+1
		Q=zeros(x)	
		#println(size(Q))
		for k in 1:nc
			theta::Float64=-1.0*(pi/2)*(gauss_ab8[k]-1)
			zk::Complex{Float64}=(lmin+lmax)/2+r*exp(0+im*theta)
			#G=zk*speye(n)-\(B,A)
			G=zk*B-A
					
			Qk=zeros(Complex128,n,m0)

			#=Gr=full(G)
			br=convert(Array{Complex128,2},x)
			Bbr=B*br
			Qk=improveSolve(Gr,Bbr,10,1)=#
		
			bx=B*x
			Qk=\(G,bx)

			println(norm(bx-G*Qk))

			#Qk=\(G,x)
		
			#for i in 1:m0	
			#Qk[:,i]=solveSys(G,x[:,i],x[:,i],1)
			#end

			#x0=rand(n,m0)
			#Qk=solveSys3(G,x,x,2,3)

			p=1	
			delta=norm(G*Qk[:,p]-B*x[:,p])/norm(B*x[:,p])
			println("        Lin sys error = $delta")
			
			Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
		end

		#(QQ,RR)=qr(Q)
		#Q[:,1:m0]=QQ[:,1:m0]
		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*(B*Q))		
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		
		resvecs=A*x-B*x*spdiagm(lest)
		reslist=zeros(m0,1)
		for i in 1:m0
			reslist[i]=norm(resvecs[:,i])/norm(x[:,i])
		end
		(res,ind)=findmin(reslist)
		
		println("FEAST Iteration $it: res = $res")
		println("~~~X norm = $(norm(x))")
		if isnan(res)
			#println("Got NaN! Printin X...")
			#println(x)
			break
		end
		residuals[it]=res
		#println(lest)
		#println(res)
	end
	return (lest,x,residuals)
end

end
