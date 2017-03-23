#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 ellipse"
#$bindir/time_sddriver $matdir/Na5 Na5_ellipse/Na5 >Na5_ellipse/Na5.out
#$bindir/plot_contour Na5_ellipse/Na5 Na5_ellipse/Na5>/dev/null
echo "Na5 line"
#$bindir/time_sddriver $matdir/Na5 Na5_line/Na5 >Na5_line/Na5.out
#$bindir/plot_contour Na5_line/Na5 Na5_line/Na5>/dev/null
echo "Na5 ellipse mid"
#$bindir/time_sddriver $matdir/Na5 Na5_ellipse_mid/Na5 >Na5_ellipse_mid/Na5.out
echo "Na5 line mid"
#$bindir/time_sddriver $matdir/Na5 Na5_line_mid/Na5 >Na5_line_mid/Na5.out

echo "Na5 ellipse 1"
#$bindir/time_sddriver $matdir/Na5 Na5_ellipse1/Na5 > Na5_ellipse1/Na5.out
#$bindir/convrate $matdir/Na5 Na5_ellipse1/Na5
#$bindir/plot_contour Na5_ellipse1/Na5 Na5_ellipse1/Na5>/dev/null
echo "Na5 ellipse 2"
#$bindir/time_sddriver $matdir/Na5 Na5_ellipse2/Na5 > Na5_ellipse2/Na5.out
#$bindir/convrate $matdir/Na5 Na5_ellipse2/Na5
#$bindir/plot_contour Na5_ellipse2/Na5 Na5_ellipse2/Na5>/dev/null
echo "Na5 ellipse 3"
#$bindir/time_sddriver $matdir/Na5 Na5_ellipse3/Na5 > Na5_ellipse3/Na5.out
#$bindir/convrate $matdir/Na5 Na5_ellipse3/Na5
#$bindir/plot_contour Na5_ellipse3/Na5 Na5_ellipse3/Na5>/dev/null
echo "Na5 ellipse 4"
$bindir/time_sddriver $matdir/Na5 Na5_ellipse4/Na5 > Na5_ellipse4/Na5.out
$bindir/convrate $matdir/Na5 Na5_ellipse4/Na5
$bindir/plot_contour Na5_ellipse4/Na5 Na5_ellipse4/Na5>/dev/null
echo "Na5 ellipse 5"
$bindir/time_sddriver $matdir/Na5 Na5_ellipse5/Na5 > Na5_ellipse5/Na5.out
$bindir/convrate $matdir/Na5 Na5_ellipse5/Na5
$bindir/plot_contour Na5_ellipse5/Na5 Na5_ellipse5/Na5>/dev/null
