set terminal pdfcairo enhanced size 5.00in,1.75in
set output "unevencol.pdf"
set xzeroaxis
set grid xtics
set grid ytics
unset key
set xlabel "Real Part"
set ylabel "Imaginary Part"
plot  "xout.dat" w p lt 1 pt 2 ps 2 lw 3,"xin.dat" w p lt 2 pt 13 ps 2 lw 3

set output "unevencolrho.pdf"
plot "newxout.dat" w p lt 1 pt 2 ps 2 lw 3,"newxin.dat" w p lt 2 pt 13 ps 2 lw 3 

