# This is a script for Gnuplot 4.0.

reset
##set terminal x11 reset
##set terminal x11 0
set terminal epslatex color
set output "dz_ds.eps"

b(s) = exp(1 - 1/s**2)
pa = 5e4
a(s) = pa * (s - b(s))
ps = 101325.
p(s) = a(s) + b(s) * ps

H = 7.
dz_ds(s) = -H * (pa + (ps - pa) * 2 / s**3 * b(s)) / p(s)

unset key
set xlabel "$s$"
set ylabel "\\shortstack{$z'(s)$ \\\\ (km)}"
plot [s=0.2:1] [-30:-10] dz_ds(s)

set output
