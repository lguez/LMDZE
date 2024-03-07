# This is a script for Gnuplot. This script plots the coefficients
# "ap" and "bp" used for the vertical grid in LMDZ.

reset
##set terminal x11 reset
set terminal epslatex color
set output "disvert.eps"

pa = 5e4
ps = 101325.
bp(s) = exp(1. - 1. / s**2)
ap(s) = pa * (s - bp(s))
p(s) = ap(s) + bp(s) * ps

##set terminal x11 0 title "bp"
set parametric
##set xlabel "bp"
set ylabel "$s$"
set yrange [] reverse
set dummy s
set trange [0.01:1]
##unset key
##plot bp(s), s

# # set terminal x11 1 title "ap"
# # set xlabel "ap"
# # plot ap(s), s

##set terminal x11 2 title "p"
##unset xlabel
##set key default
plot p(s), s title "$p$", ap(s), s title "$a$", \
     bp(s) * ps, s title "$b \\times p_s$"

set output
