#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 flat mid"
$bindir/time_sddriver $matdir/Na5 Na51_flat_mid/Na5 >Na51_flat_mid/Na5.out
$bindir/plot_contour Na51_flat/Na5 Na51_flat_mid/Na5 >/dev/null

