set terminal svg background rgb "#ffffff"
set grid

set output "scaling_time.svg"
set title "Jacobi diagonalization timing"
set xlabel "Matrix size N"
set ylabel "Time [s]"
plot "scaling.times.txt" using 1:2 with linespoints lw 2 pt 7 title "measured"

set output "scaling_ratio.svg"
set title "Jacobi scaling check: t(N)/N^3"
set xlabel "Matrix size N"
set ylabel "t/N^3 [s]"
plot "scaling.ratio.txt" using 1:3 with linespoints lw 2 pt 7 title "t/N^3"
