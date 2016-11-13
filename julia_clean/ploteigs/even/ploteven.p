set terminal pdfcairo enhanced size 5.00in,1.75in
set output "even.pdf"
set xzeroaxis
set grid xtics
set grid ytics
unset key
set xlabel "Real Part"
set ylabel "Imaginary Part"
plot "xin.dat" w p pt 2 ps 2 lw 3, "xout.dat" w p lt 1 pt 2 ps 2 lw 3
