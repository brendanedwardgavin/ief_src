

module juliaFEAST

export feast,feastG

include("gmresRP.jl")
include("lsSolve.jl")

gauss_weight8=[0.362683783378361;0.362683783378361;0.313706645877887;0.313706645877887;0.222381034453374;0.222381034453374;0.101228536290376;0.101228536290376]
gauss_ab8=[0.183434642495649;-0.183434642495649;0.525532409916328;-0.525532409916328;0.796666477413626;-0.796666477413626;0.960289856497536;-0.960289856497536]

function feastfilt(A,Y,lmin,lmax)
	nc=8
	(n,m)=size(Y)

	r=(lmax-lmin)/2.0
	
	
	Q=zeros(n,m)	
	
	for k in 1:nc
		theta=-1.0*(pi/2)*(gauss_ab8[k]-1)
		zk=(lmin+lmax)/2+r*exp(0+im*theta)
		G=zk*speye(n)-A
		Qk=\(G,Y)
		Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
	end

	return Q	
end

function feast(A::SparseMatrixCSC{Float64,Int64},x0::Array{Float64,2},lmin::Float64,lmax::Float64,m0::Int64,eps::Float64,maxit::Int64)
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
			G=zk*speye(n)-A
			

			#Qk=\(G,x)
			#=if it>1
				alph=-0.00
				#println(lest)
				lzk=diagm(1./((zk-lest)))
				lzk2=diagm(zk-lest)

				Qk=x #(initial guess)
				for i in 1:maxlinit
					Qk=alph*Qk+(1-alph)*(x*lzk-x*(lzk*(x'*(G*Qk)))-x*(x'*Qk))
				end

				#Qk=\(G,x)
			else
				Qk=\(G,x)
			end=#

			
			Qk=zeros(Complex128,n,m0)
			#=for i in 1:m0
				Gr=[real(G) -imag(G);imag(G) real(G)]
				b=x[:,i]
				br=[real(x[:,i]);zeros(Float64,n)]
				x0=zeros(Float64,2*n)
				x0[:]=br

				M=eye(Float64,2*n)
				
				#=if it>1
				lzk2=diagm(zk-lest)
				lzk=diagm(1./((zk-lest)))
				M=x*lzk*x'

				M=[real(M) -imag(M);imag(M) real(M)]
				end=#
		
				restrt=5
				gmmaxit=1
				tol=1e-16
				
				x0=gmresRP(Gr,x0,br,M,restrt,gmmaxit,tol)
				
				#sol=solveSys(Gr,br,br,10)
				
				Qk[:,i]=x0[1:n]+im*x0[n+1:2*n]
				#Qk[:,i]=sol[1:n]+im*sol[n+1:2*n]
			end=#

			Gr=full(G)
			br=convert(Array{Complex128,2},x)
			Qk=improveSolve(Gr,br,10,1)


			#Qk=\(G,x)
		
			#for i in 1:m0	
			#Qk[:,i]=solveSys(G,x[:,i],x[:,i],1)
			#end

			#x0=rand(n,m0)
			#Qk=solveSys3(G,x,x,2,3)

			p=1	
			delta=norm(G*Qk[:,p]-x[:,p])/norm(x[:,p])
			println("        Lin sys error = $delta")
			
			Q=Q-(gauss_weight8[k]/2)*real(r*exp(im*theta)*Qk)
		end

		#(QQ,RR)=qr(Q)
		#Q[:,1:m0]=QQ[:,1:m0]
		for i in 1:m0
			Q[:,i]=Q[:,i]/norm(Q[:,i])
		end

		Aq=full(Q'*(A*Q))
		Bq=full(Q'*Q)		
		F=eigfact(Aq,Bq)
		lest=F[:values]
		xest=F[:vectors]
		x=Q*xest

		
		resvecs=A*x-x*spdiagm(lest)
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
			G=zk*B-A
			
			#Qk=\(G,x)
			if 1==1
			Qk=zeros(Complex128,n,m0)
			Gr=full(G)
			#Gr=full(\(B,Gr))
			br=convert(Array{Complex128,2},x)
			Qk=improveSolveG(Gr,full(B),br,20,1)
			end

			#Qk=\(G,B*x)
		
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
			reslist[i]=norm(resvecs[:,i])/norm(B*x[:,i])
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
