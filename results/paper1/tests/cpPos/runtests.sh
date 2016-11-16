#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 ellipse"
$bindir/time_sddriver $matdir/Na5 Na5_ellipse/Na5 >Na5_ellipse/Na5.out
echo "Na5 line"
$bindir/time_sddriver $matdir/Na5 Na5_line/Na5 >Na5_line/Na5.out

