flags=-fpp -mkl -traceback #-check all
all:	
	ifort ../feast_gmres.f90 test.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a $(flags) -o ./test
