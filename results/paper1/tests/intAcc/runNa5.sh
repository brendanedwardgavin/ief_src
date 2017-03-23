#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 gmres 16"
$bindir/time_sddriver $matdir/Na5 Na5_gm16/Na5 >Na5_gm16/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm16/Na5
$bindir/plot_contour Na5_gm16/Na5 Na5_gm16/Na5
echo "Na5 gmres 4"
$bindir/time_sddriver $matdir/Na5 Na5_gm4/Na5 >Na5_gm4/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm4/Na5
$bindir/plot_contour Na5_gm4/Na5 Na5_gm4/Na5

echo "Na5 gmres 8"
$bindir/time_sddriver $matdir/Na5 Na5_gm8/Na5 >Na5_gm8/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm8/Na5  
$bindir/plot_contour Na5_gm8/Na5 Na5_gm8/Na5



echo "Na5 gmres 16 flat"
$bindir/time_sddriver $matdir/Na5 Na5_gm16_flat/Na5 >Na5_gm16_flat/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm16_flat/Na5
$bindir/plot_contour Na5_gm16_flat/Na5 Na5_gm16_flat/Na5

echo "Na5 gmres 4 flat"
$bindir/time_sddriver $matdir/Na5 Na5_gm4_flat/Na5 >Na5_gm4_flat/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm4_flat/Na5
$bindir/plot_contour Na5_gm4_flat/Na5 Na5_gm4_flat/Na5

echo "Na5 gmres 8 flat"
$bindir/time_sddriver $matdir/Na5 Na5_gm8_flat/Na5 >Na5_gm8_flat/Na5.out
$bindir/convrate $matdir/Na5 Na5_gm8_flat/Na5  
$bindir/plot_contour Na5_gm8_flat/Na5 Na5_gm8_flat/Na5



echo "Na5 pardiso 16"
#$bindir/time_sddriver $matdir/Na5 Na5_pa16/Na5 >Na5_pa16/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa16/Na5  
echo "Na5 pardiso 4"
#$bindir/time_sddriver $matdir/Na5 Na5_pa4/Na5 >Na5_pa4/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa4/Na5 
echo "Na5 pardiso 8"
#$bindir/time_sddriver $matdir/Na5 Na5_pa8/Na5 >Na5_pa8/Na5.out
#$bindir/convrate $matdir/Na5 Na5_pa8/Na5 
