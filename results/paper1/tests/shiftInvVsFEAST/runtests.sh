#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 flat"
#$bindir/time_sddriver $matdir/Na5 Na51_flat/Na5 >Na51_flat/Na5.out
#$bindir/plot_contour Na51_flat/Na5 Na51_flat/Na5 >/dev/null
echo "Na5 ellipse"
#$bindir/time_sddriver $matdir/Na5 Na51_ellipse/Na5 >Na51_ellipse/Na5.out
#$bindir/plot_contour Na51_ellipse/Na5 Na51_ellipse/Na5 >/dev/null

echo "Na5 flat mid"
$bindir/time_sddriver $matdir/Na5 Na51_flat_mid/Na5 >Na51_flat_mid/Na5.out
$bindir/plot_contour Na51_flat/Na5 Na51_flat_mid/Na5 >/dev/null

echo "Na5 ellipse mid"
#$bindir/time_sddriver $matdir/Na5 Na51_ellipse_mid/Na5 >Na51_ellipse_mid/Na5.out
#$bindir/plot_contour Na51_ellipse_mid/Na5 Na51_ellipse_mid/Na5 >/dev/null

echo "Na5 above mid"
$bindir/time_sddriver $matdir/Na5 Na51_above_mid/Na5 >Na51_above_mid/Na5.out
$bindir/plot_contour Na51_above_mid/Na5 Na51_above_mid/Na5 >/dev/null
