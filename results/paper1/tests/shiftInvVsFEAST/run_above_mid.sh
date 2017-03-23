#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 above mid"
$bindir/time_sddriver $matdir/Na5 Na51_above_mid/Na5 >Na51_above_mid/Na5.out
$bindir/plot_contour Na51_above_mid/Na5 Na51_above_mid/Na5 >/dev/null
