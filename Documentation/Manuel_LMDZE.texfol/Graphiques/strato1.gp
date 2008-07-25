# This is a script for Gnuplot 4.0.

reset
##set term x11 reset; set term x11 0
set terminal epslatex color; set output "strato1.eps"

f(x, y) = tanh(x) + 0.5 * tanh((x - y) / 2)
llm = 70
set xlabel "$x$"
set ylabel "$f(x;\\nombre{11,6})$"
unset key
plot [-11./5.: (llm - 12.) / 5.] f(x, (llm - 12.) / 5.)

set output
