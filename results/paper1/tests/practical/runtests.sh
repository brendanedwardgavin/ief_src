#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "CO gmres"
#$bindir/time_sddriver $matdir/CO CO_gm/CO > CO_gm/CO.out
echo "CO pardiso"
#$bindir/time_sddriver $matdir/CO CO_pa/CO > CO_pa/CO.out
echo "Ga3As3H12 gmres"
$bindir/time_sddriver $matdir/Ga3As3H12 Ga3As3H12_gm/Ga3As3H12 > Ga3As3H12_gm/Ga3As3H12.out
echo "Ga3As3H12 pardiso"
$bindir/time_sddriver $matdir/Ga3As3H12 Ga3As3H12_pa/Ga3As3H12 > Ga3As3H12_pa/Ga3As3H12.out
echo "Na5 gmres"
#$bindir/time_sddriver $matdir/Na5 Na5_gm/Na5 > Na5_gm/Na5.out
echo "Na5 Pardiso"
#$bindir/time_sddriver $matdir/Na5 Na5_pa/Na5 > Na5_pa/Na5.out
echo "Qe4213 gmres"
#$bindir/time_sddriver $matdir/H_matrix_4213 qe4213_gm/H_matrix_4213 > qe4213_gm/H_matrix_4213.out
#echo "Qe4213 pardiso"
#$bindir/time_sddriver $matdir/H_matrix_4213 qe4213_pa/H_matrix_4213 > qe4213_pa/H_matrix_4213.out

