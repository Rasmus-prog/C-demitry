set terminal svg background rgb "#ffffff"
set output "convergence.svg"

set multiplot layout 1,2 title "Hydrogen s-wave convergence of epsilon_0"

set key top right
set grid

set xlabel "dr"
set ylabel "epsilon_0 (Hartree)"
set title "Fixed rmax: epsilon_0 vs dr"
set yrange [-0.505:-0.47]
set xrange [*:*] reverse
plot \
  "convergence_dr.dat" using 1:2 with linespoints lw 2 pt 7 title "numeric", \
  "convergence_dr.dat" using 1:3 with lines lw 2 title "exact"
set xrange [*:*] noreverse

set xlabel "rmax"
set ylabel "epsilon_0 (Hartree)"
set title "Fixed dr: epsilon_0 vs rmax"
set yrange [-0.505:-0.46]
plot \
  "convergence_rmax.dat" using 1:2 with linespoints lw 2 pt 7 title "numeric", \
  "convergence_rmax.dat" using 1:3 with lines lw 2 title "exact"

unset multiplot