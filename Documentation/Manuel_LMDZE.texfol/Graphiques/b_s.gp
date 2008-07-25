# This is a script for Gnuplot 4.0.

reset
##set term x11 reset
##set term x11 0
set terminal epslatex color
set output "b_s.eps"

b(s) = exp(1 - 1 / s**2)
set logscale  y
set xlab "s"
set ylab "$b(s)/s$"
unset key
plot [s=0.3:1] b(s) / s

set output
