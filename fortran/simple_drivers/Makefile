all:	
	ifort ../feast_gmres.f90 zfeastgmres_driver.f90 ../../myFeast3.0/lib/x64/libfeast.a -mkl -traceback -check -gen-interfaces -warn -fpe0 -warn nounused -o zfeastgmres_driver
	ifort ../feast_gmres.f90 dfeastgmres_driver.f90 ../../myFeast3.0/lib/x64/libfeast.a -mkl -traceback -check -gen-interfaces -warn -fpe0 -warn nounused -o dfeastgmres_driver

