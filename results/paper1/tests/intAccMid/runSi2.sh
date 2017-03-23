#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Si2 minres 2"
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr2/Si2 >Si2_mid_mr2/Si2.out
#$bindir/convrate $matdir/Na5 Na5_mr2/Na5  
#$bindir/plot_contour Na5_mr2/Na5 Na5_mr2/Na5

echo "Na5 gmres 2"
#$bindir/time_sddriver $matdir/Na5 Na5_gm2/Na5 >Na5_gm2/Na5.out
#$bindir/convrate $matdir/Na5 Na5_gm2/Na5  
#$bindir/plot_contour Na5_gm2/Na5 Na5_gm2/Na5



echo "Na5 pardiso 8"
#$bindir/time_sddriver $matdir/Na5 Na5_pa8/Na5 >Na5_pa8/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa8/Na5  
echo "Na5 pardiso 4"
#$bindir/time_sddriver $matdir/Na5 Na5_pa4/Na5 >Na5_pa4/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa4/Na5 
echo "Na5 pardiso 2"
#$bindir/time_sddriver $matdir/Na5 Na5_pa2/Na5 >Na5_pa2/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa2/Na5 
