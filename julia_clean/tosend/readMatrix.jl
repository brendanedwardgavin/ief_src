function dreadMatrixCoo(filename)

	file=open(filename)
	
	#skip comments
	line=readline(file)
	while length(chomp(line))==0 || (length(line)>0 && line[1]=='%')
		line=readline(file)
	end

	flds=split(line)

	rows=parse(Int64,flds[1])
	cols=parse(Int64,flds[2])
	nnz=parse(Int64,flds[3])

	rc=Array(Int64, nnz)
	cc=Array(Int64, nnz)
	anz=Array(Float64, nnz)
	
	for i in 1:nnz

		flds=split(readline(file))	
		rc[i]=parse(Int64,flds[1])
		cc[i]=parse(Int64,flds[2])
		anz[i]=parse(Float64,flds[3])

	end

	mat=sparse(rc,cc,anz)

	return mat
end
