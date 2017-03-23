#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Si2 mid 2"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr2/Si2 >Si2_mid_mr2/Si2.out
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr2_20/Si2 >Si2_mid_mr2_20/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_mr2/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_mr2/Si2  
#$bindir/plot_contour Si2_mid_mr2/Si2 Si2_mid_mr2/Si2 >/dev/null

echo "Si2 mid 4"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2.out
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr4_20/Si2 >Si2_mid_mr4_20/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_mr4/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_mr4/Si2  
#$bindir/plot_contour Si2_mid_mr4/Si2 Si2_mid_mr4/Si2 >/dev/null

echo "Si2 mid 8"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_mr8/Si2 >Si2_mid_mr8/Si2.out
$bindir/time_sddriver $matdir/Si2 Si2_mid_mr8_20/Si2 >Si2_mid_mr8_20/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_mr8/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_mr8/Si2  
#$bindir/plot_contour Si2_mid_mr8/Si2 Si2_mid_mr8/Si2 >/dev/null


echo "Si2 par 2"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_pa2/Si2 >Si2_mid_pa2/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_pa2/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_pa2/Si2  
#$bindir/plot_contour Si2_mid_pa2/Si2 Si2_mid_pa2/Si2 >/dev/null

echo "Si2 mid 4"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_pa4/Si2 >Si2_mid_pa4/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_pa4/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_pa4/Si2  
#$bindir/plot_contour Si2_mid_pa4/Si2 Si2_mid_pa4/Si2 >/dev/null

echo "Si2 mid 8"
#$bindir/time_sddriver $matdir/Si2 Si2_mid_pa8/Si2 >Si2_mid_pa8/Si2.out
#$bindir/plot_alphas_empirical $matdir/Si2 Si2_mid_mr4/Si2 >Si2_mid_mr4/Si2alphas.out
#$bindir/plot_convrate $matdir/Si2 Si2_mid_pa8/Si2 
#$bindir/convrate $matdir/Si2 Si2_mid_pa8/Si2  
#$bindir/plot_contour Si2_mid_pa8/Si2 Si2_mid_pa8/Si2 >/dev/null

