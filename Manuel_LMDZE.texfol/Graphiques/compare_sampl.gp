# This is a script for Gnuplot 4.0.

reset
set term x11 reset
set term x11 0

set style data linespoints
set yrange [] reverse
set logscale y
set mytics 10
set ylabel "\\shortstack{pression \\\\ (Pa)}"
set xlabel "indice"
set grid ytics
##set pointsize 2
set key left

plot "test_disvert_LMD5.csv" using 4 title "LMD5", "test_disvert_param.csv" using 4 title "param", "test_disvert_strato1.csv" using 4 title "strato1", "test_disvert_strato2.csv" using 4 title "strato2"