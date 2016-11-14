include("geneigs.jl")

filename=ARGS[1]

eigdist(x)=1.0 #uniform distribution

emin=-3.0
emax=5.0

n=1000

println("Generating eigenvalues")
eigvals=geneigs(n,eigdist,emin,emax)

println("Doing QR")
X=rand(n,n)
(Q,R)=qr(X)

println("forming matrix")
A=Q*diagm(eigvals)*Q'

println("writing matrix")
open(filename,"w") do f
    write(f,"$n $n $(n*n)\n")
    for i in 1:n
    if (mod(i,100)==0)
        println("row $i")
    end
    
    for j in 1:n
        write(f,"$i $j $(A[i,j])\n")
    end
    end
end
