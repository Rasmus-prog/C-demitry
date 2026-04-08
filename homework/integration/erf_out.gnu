set terminal svg background rgb "white"
set output "erf_out.svg"


set title "erf(x)"
set logscale x
set xlabel "acc"
set ylabel "erf(x)"
set title "erf(x) as a function of acc"
plot "erf_out.data" with linespoints title "erf(x)"