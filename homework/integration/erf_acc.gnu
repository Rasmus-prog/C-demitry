set terminal svg background rgb "white"
set output "erf_acc.svg"


set title "abs(erf(1) - exact) vs acc"
set logscale x
set xlabel "acc"
set ylabel "abs(erf(1) - exact)"
set title "erf(1) as a function of acc"
plot "erf_acc.data" with linespoints title "abs(erf(1) - exact) vs acc"