# This is a script for Gnuplot. This script plots the coefficients
# "ap" and "bp" used for the vertical grid in LMDZ.

reset
set terminal x11 reset

pa = 5e4
ps = 101325.
bp(sigma) = exp(1. - 1. / sigma**2)
ap(sigma) = pa * (sigma - bp(sigma))
p(sigma) = ap(sigma) + bp(sigma) * ps

set terminal x11 0 title "bp"
set parametric
set xlabel "bp"
set ylabel "sigma"
set yrange [] reverse
set dummy sigma
set trange [0.01:1]
unset key
plot bp(sigma), sigma

set terminal x11 1 title "ap"
set xlabel "ap"
plot ap(sigma), sigma

set terminal x11 2 title "p"
unset xlabel
set key default
plot p(sigma), sigma title "p", ap(sigma), sigma title "ap", bp(sigma) * ps, \
     sigma title "bp * ps"
