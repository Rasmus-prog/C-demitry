# set terminal pngcairo size 1200,500
# set output 'convergence.svg'
set terminal svg background rgb "#ffffff"
set output "convergence.svg"

set multiplot layout 1,2 title 'Hydrogen s-wave convergence of epsilon_0'

set grid
set key top right

set xlabel 'dr (Bohr)'
set ylabel 'epsilon_0 (Hartree)'
set title 'Convergence vs dr (fixed rmax)'
set xrange [0.5:0.1]
set yrange [-0.505:-0.47]
plot 'convergence_dr.dat' using 1:2 with linespoints lw 2 pt 7 title 'numeric epsilon_0', \
     -0.5 with lines lw 2 dt 2 title 'exact -1/2'

set xlabel 'rmax (Bohr)'
set ylabel 'epsilon_0 (Hartree)'
set title 'Convergence vs rmax (fixed dr)'
set xrange [*:*]
set yrange [-0.505:-0.47]
plot 'convergence_rmax.dat' using 1:2 with linespoints lw 2 pt 7 title 'numeric epsilon_0', \
     -0.5 with lines lw 2 dt 2 title 'exact -1/2'

unset multiplot
unset output
