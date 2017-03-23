#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 low 75"
#$bindir/time_sddriver $matdir/Na5 Na5low_75/Na5 >Na5low_75/Na5.out
echo "Na5 low 100"
#$bindir/time_sddriver $matdir/Na5 Na5low_100/Na5 >Na5low_100/Na5.out
echo "Na5 low 200"
#$bindir/time_sddriver $matdir/Na5 Na5low_200/Na5 >Na5low_200/Na5.out
echo "Na5 low 800"
#$bindir/time_sddriver $matdir/Na5 Na5low_800/Na5 >Na5low_800/Na5.out
echo "Na5 low 787"
$bindir/time_sddriver $matdir/Na5 Na5low_787/Na5 >Na5low_787/Na5.out


echo "Arnoldi low"
#$bindir/svd_arpack $matdir/Na5 L SA 50 800 >arnoldilow/Na5_800.out
#$bindir/svd_arpack $matdir/Na5 L SA 50 75 >arnoldilow/Na5_75.out
#$bindir/svd_arpack $matdir/Na5 L SA 50 100 >arnoldilow/Na5_100.out
#$bindir/svd_arpack $matdir/Na5 L SA 50 200 >arnoldilow/Na5_200.out
#$bindir/svd_arpack $matdir/Na5 L SA 50 787 >arnoldilow/Na5_787.out
