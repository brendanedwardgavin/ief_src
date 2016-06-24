flags=-fpp -mkl -traceback -check all
all:	
	ifort ../feast_gmres.f90 time_szfeastgmres.f90 time_szdriver.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a $(flags)  -o ../../bin/time_szdriver
	ifort ../feast_gmres.f90 time_szfeastgmres.f90 time_sddriver.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a $(flags) -o ../../bin/time_sddriver
	ifort ../feast_gmres.f90 time_dzfeastgmres.f90 time_dddriver.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a $(flags) -o ../../bin/time_dddriver
	ifort lapack_dddriver.f90 $(flags) -o ../../bin/lapack_dddriver
