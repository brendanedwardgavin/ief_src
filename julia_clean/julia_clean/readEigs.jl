function readEigs(filename)

	srand(2)
	file=open(filename)
	a=readline(file)
	n=parse(Int64,a)
	eigs=zeros(n)
	for i in 1:n
		a=readline(file)
		eigs[i]=parse(Float64,a)
	end

	x=rand(n,n)
	(q,r)=qr(x)
	x=q*diagm(eigs)*q'

	return x #real symmetric matrix with eigs eigenvalues
end

function zreadMat(filename)
        file=open(filename)
        a=readline(file)
        n=parse(Int64,a)
        A=zeros(Complex128,n,n)
        for i in 1:n
	for j in 1:n
                a=readline(file)
		z=readdlm(IOBuffer(a))
                A[i,j]=z[1]+z[2]*im
	end
	num=rand()
	if num>0.95
		println("$(round(100*i/n,1))%")
	end
        end

	#println(A[2,1])

        return A
end
