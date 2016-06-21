all:	
	ifort ../feast_gmres.f90 time_szfeastgmres.f90 time_szdriver.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a -fpp -mkl -traceback -o time_szdriver
	ifort ../feast_gmres.f90 time_szfeastgmres.f90 time_sddriver.f90 ../../myFeast3.0/lib/x64/libfeast.a ../../myFeast3.0/lib/x64/libfeast_sparse.a -fpp -mkl -traceback -o time_sddriver
