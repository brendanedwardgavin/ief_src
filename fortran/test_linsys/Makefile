flags=-mkl -traceback -g -check all
feastdir=/opt/FEAST/3.0/lib/x64
bindir=../../bin

all:	
	ifort ../rundata.f90 ../zminres.f90 ../feast_gmres.f90 dlinsys.f90 $(feastdir)/libfeast.a $(flags)  -o $(bindir)/dlinsys	
