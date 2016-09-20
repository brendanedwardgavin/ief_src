#plot "polysamples.dat" w l lw 2, "polyzeros.dat" w p ps 2 lw 3
plot "truepolyf.dat" w l lw 2, "kpolyf.dat" w l lw 2, "trueeigs.dat" u 1:2 w p ps 2 lw 3, "" u 1:2:3 w labels center offset 0,1 , "approxeigs.dat" w p ps 2 lw 3
