#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Si2 mid 20"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr20/Si2 >Si2_mid_mr20/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr20/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr20/Si2  
$bindir/plot_contour Si2_mid_mr20/Si2 Si2_mid_mr20/Si2 >/dev/null

echo "Si2 mid 21"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr21/Si2 >Si2_mid_mr21/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr21/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr21/Si2  
$bindir/plot_contour Si2_mid_mr21/Si2 Si2_mid_mr21/Si2 >/dev/null

echo "Si2 mid 30"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr30/Si2 >Si2_mid_mr30/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr30/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr30/Si2  
$bindir/plot_contour Si2_mid_mr30/Si2 Si2_mid_mr30/Si2 >/dev/null

echo "Si2 mid 40"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr40/Si2 >Si2_mid_mr40/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr40/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr40/Si2  
$bindir/plot_contour Si2_mid_mr40/Si2 Si2_mid_mr40/Si2 >/dev/null

echo "Si2 mid 200"
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr200/Si2 >Si2_mid_mr200/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
$bindir/plot_convrate $matdir/Si2 Si2_mid_mr40/Si2 
$bindir/convrate $matdir/Si2 Si2_mid_mr40/Si2  
$bindir/plot_contour Si2_mid_mr40/Si2 Si2_mid_mr40/Si2 >/dev/null
