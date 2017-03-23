#!/bin/bash
bindir=../../../../bin
matdir=../../matrices

echo "Na5 pardiso"
#$bindir/time_sddriver $matdir/Na5 Na5Pardiso/Na5 >Na5Pardiso/Na5.out
echo "Na5 1"
#$bindir/time_sddriver $matdir/Na5 Na51/Na5 >Na51/Na5.out
echo "Na5 2"
#$bindir/time_sddriver $matdir/Na5 Na52/Na5 >Na52/Na5.out
echo "Na5 3"
#$bindir/time_sddriver $matdir/Na5 Na53/Na5 >Na53/Na5.out

echo "Na5 pardiso mid"
#$bindir/time_sddriver $matdir/Na5 Na5Pardiso_mid/Na5 >Na5Pardiso_mid/Na5.out
echo "Na5 1 mid"
#$bindir/time_sddriver $matdir/Na5 Na51_mid/Na5 >Na51_mid/Na5.out
echo "Na5 2 mid"
#$bindir/time_sddriver $matdir/Na5 Na52_mid/Na5 >Na52_mid/Na5.out
echo "Na5 3 mid"
#$bindir/time_sddriver $matdir/Na5 Na53_mid/Na5 >Na53_mid/Na5.out

echo "Na5 pardiso 8  mid"
$bindir/time_sddriver $matdir/Na5 cp8/Na5Pardiso_mid/Na5 >cp8/Na5Pardiso_mid/Na5.out
echo "Na5 1 8  mid"
$bindir/time_sddriver $matdir/Na5 cp8/Na51_mid/Na5 >cp8/Na51_mid/Na5.out
echo "Na5 2 8  mid"
$bindir/time_sddriver $matdir/Na5 cp8/Na52_mid/Na5 >cp8/Na52_mid/Na5.out
echo "Na5 3 8  mid"
$bindir/time_sddriver $matdir/Na5 cp8/Na53_mid/Na5 >cp8/Na53_mid/Na5.out
