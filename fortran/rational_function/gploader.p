set multiplot layout 1,2
set xrange [-1.0:1.0]
set yrange [1e-6:1.5]
#set logscale y
set format y "%12.1f"
unset key
set title "Rational Function" font "Arial,18"
set xlabel "{/Symbol l}" font "Arial,15"
set ylabel "{/Symbol r}({/Symbol l})" offset 9 font "Arial,15"
plot sprintf("output/%i.dat",i) u 1:2 w l lw 2
unset logscale

set key
set title sprintf("Contour, Ellipse = %i %%",i) font "Arial,18" 
set xrange [-1.0:1.0]
set yrange [-0.75:0.75]
set zeroaxis lw 2
set ylabel "Im(z)"
set xlabel "Re(z)"
plot sprintf("output/%iel.dat",i) w l lt 7 lw 2 title "Contour",sprintf("output/%icp.dat",i) w p ps 2 pt 13 lc 2 title "Quadrature Points"
unset zeroaxis
unset multiplot
i=i+1
if(i<n) reread
