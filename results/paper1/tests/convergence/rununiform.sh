#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Uniform pardiso"
#$bindir/time_sddriver $matdir/uniform uniformPardiso/uniform>uniformPardiso/uniform.out
echo "Uniform 1"
$bindir/time_sddriver $matdir/uniform uniform1/uniform>uniform1/uniform.out
echo "Uniform 2"
$bindir/time_sddriver $matdir/uniform uniform2/uniform>uniform2/uniform.out
echo "Uniform 3"
$bindir/time_sddriver $matdir/uniform uniform3/uniform>uniform3/uniform.out

