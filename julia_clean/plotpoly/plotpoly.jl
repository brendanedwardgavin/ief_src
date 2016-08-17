include("../geneigs.jl")
include("writefiles.jl")

#eigdist(x)=exp(-(x-10.0)^2)+exp(-(x+10.0)^2)
eigdist(x)=1.0
emin=-12.0
emax=12.0
n=10

trueeigs=geneigs(n,eigdist,emin,emax)

for i in 1:n
	println("$i $(trueeigs[i])")
end

function poly(x)
	fx=1.0
	for i in 1:n
		fx=fx*(x-trueeigs[i])	
	end

	return fx
end

de=abs(emax-emin)

writefunction(poly,500,minimum(trueeigs),maximum(trueeigs))

writezeros(trueeigs)
