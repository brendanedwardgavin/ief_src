#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Uniform 1"
#$bindir/time_sddriver $matdir/uniform uniform1/uniform>uniform1/uniform.out
echo "Uniform 2"
#$bindir/time_sddriver $matdir/uniform uniform2/uniform>uniform2/uniform.out
echo "Uniform 3"
#$bindir/time_sddriver $matdir/uniform uniform3/uniform>uniform3/uniform.out
echo "Uniform 4"
$bindir/time_sddriver $matdir/uniform uniform4/uniform>uniform4/uniform.out
echo "Uniform 5"
$bindir/time_sddriver $matdir/uniform uniform5/uniform>uniform5/uniform.out


