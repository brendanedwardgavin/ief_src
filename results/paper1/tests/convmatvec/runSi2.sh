#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Si2 mid 4"
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr4/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr4/Si2  
$bindir/plot_contour Si2_mid_mr4/Si2 Si2_mid_mr4/Si2 >/dev/null


echo "Si2 low 4"
#$bindir/time_sddriver $matdir/Si2 Si2_low_mr4/Si2 >Si2_low_mr4/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_low_mr4/Si2 >Si2_low_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_low_mr4/Si2
$bindir/convrate $matdir/Si2 Si2_low_mr4/Si2  
$bindir/plot_contour Si2_low_mr4/Si2 Si2_low_mr4/Si2

echo "Si2 midtall 4"
#$bindir/time_sddriver $matdir/Si2 Si2_midtall_mr4/Si2 >Si2_midtall_mr4/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_midtall_mr4/Si2 
$bindir/convrate $matdir/Si2 Si2_midtall_mr4/Si2  
$bindir/plot_contour Si2_midtall_mr4/Si2 Si2_midtall_mr4/Si2 >/dev/null


