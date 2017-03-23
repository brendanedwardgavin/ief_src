flags=-fpp -mkl -traceback #-check all
bindir=../../bin
feastdir=/share/apps/FEAST/3.0/lib/x64

all:	
	ifort rundata.f90 ../feast_gmres.f90 ../dzlsprim.f90 ../iterative_feast.f90 time_sddriver.f90 ../../myFeast3.0/lib/x64/libfeast.a  $(flags) -o $(bindir)/time_sddriver
	ifort lapack_dddriver.f90 ../../myFeast3.0/lib/x64/libfeast.a  $(flags) -o $(bindir)/lapack_dddriver
	ifort convrate.f90 ../../myFeast3.0/lib/x64/libfeast.a  $(flags) -o $(bindir)/convrate
	ifort plot_contour.f90 $(feastdir)/libfeast_dense.a $(feastdir)/libfeast.a $(feastdir)/libfeast_sparse.a  -fpp -mkl -traceback -o $(bindir)/plot_contour
	ifort plot_convrate.f90 ../../myFeast3.0/lib/x64/libfeast.a  $(flags) -o $(bindir)/plot_convrate
	ifort gencp.f90 ../../myFeast3.0/lib/x64/libfeast.a  $(flags) -o $(bindir)/gencp
