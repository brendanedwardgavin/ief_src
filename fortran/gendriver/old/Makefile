flags=-fpp -O3 -heap-arrays -mkl  #-check all
#flags= -heap-arrays -g -fpp -mkl -traceback  #-check all
bindir=../../bin
feastdir=/opt/FEAST/3.0/lib/x64
arpackdir=../arpack/ARPACK

all:	
	ifort ../paper1/rundata.f90 iterative_feast.f90 driver_feast_sparse.f90 $(feastdir)/libfeast.a  $(flags) -o $(bindir)/svd_feast
