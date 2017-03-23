reset
set term gif enhanced animate delay 10 size 1200,600
set output "animate.gif"
n=100
#set xrange [-1.0:1.0]
#set yrange [0.0:1.3]
i=0
load "gploader.p"
set output
