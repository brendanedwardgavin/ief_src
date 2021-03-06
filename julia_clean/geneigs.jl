using Roots #for root finding

function geneigs(n,eigdist,emin,emax)

	#eigdist should be a function with l1 norm = 1

	eigs=Array{Float64}(n)
	samples=zeros(Float64,n)

	samples[1]=0.0+1.0/(n+1)
	for i in 2:n
		samples[i]=samples[i-1]+1.0/(n+1)
	end

	intacc=1e-3

	#get norm of distribution to normalize it:
	(nrm,err)=quadgk(eigdist,emin,emax,reltol=intacc,maxevals=1e4)

	#calculate and normalize cdf:
	function cdf(x)
		if emin<=x
			(out,err)=quadgk(eigdist,emin,x,reltol=intacc,maxevals=1e4)
			return out/nrm	
		else
			return 0
		end
	end

	for i in 1:n
		fmin(x)=cdf(x)-samples[i]
		eigs[i]=fzero(fmin,emin,emax)
	end

	return eigs
end
