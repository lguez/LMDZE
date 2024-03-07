# This is a script for Gnuplot 4.0.

reset
set terminal epslatex color
set output "compare_sampl_$0.eps"
##set term x11 reset persist
##set term x11 0

set style data linespoints
set yrange [] reverse
set logscale y
set mytics 10
set ylabel "\\shortstack{pression \\\\ (Pa)}"
set xlabel "indice"
set grid ytics
##set pointsize 2
set key left

plot "test_disvert_$0_LMD5.csv" using 4 title "LMD5", "test_disvert_$0_param.csv" using 4 title "param", "test_disvert_$0_strato.csv" using 4 title "strato"

set output
