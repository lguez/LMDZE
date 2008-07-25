# This is a script for Gnuplot 4.0.

reset
##set term x11 reset; set term x11 0
set terminal epslatex color; set output "strato2.eps"

g(x) = 0.5 * (1. - tanh(x / pi * 2 - 1.))
f(x) = (1. + 7. * sin(x)**2) * g(x)**2
set xlabel "$x$"
set ylabel "$f(x)$"
unset key
plot [0: pi] f(x)

set output
