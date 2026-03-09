set terminal svg background rgb "#ffffff"
set output "timing.svg"

set title "Jacobi diagonalization timing on random symmetric NxN matrices"
set xlabel "N"
set ylabel "time (s)"
set key top left
set grid
set xrange [*:*]



f(x)=a*x**3
fit f(x) "out.times.data" using 1:2 via a

plot \
  "out.times.data" using 1:2 with points pt 7 title "measured", \
  f(x) with lines lw 2 title sprintf("fit a*N^3, a=%.3e", a)
