function writezeros(zerosarray,filename)
	m=size(vec(zerosarray),1)

	open(filename,"w") do f
		for i in 1:m
			write(f,"$(zerosarray[i])  0.0  $i\n")
		end
	end
end



function writefunction(f,samples,xmin,xmax,filename)

	fx=zeros(samples+1,2)
	for i in 1:samples+1
		fx[i,1]=xmin+(i-1)*abs(xmax-xmin)/samples
		fx[i,2]=f(fx[i,1])	
	end

	fx[:,2]=fx[:,2]/maximum(abs(fx[:,2])) #normalize

	open(filename,"w") do file
		for i in 1:samples+1
			write(file,"$(fx[i,1]) $(fx[i,2])\n")
		end
	end
end
