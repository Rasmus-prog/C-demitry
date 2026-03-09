set terminal svg background rgb "#ffffff"
set output "wavefunctions.svg"

set title "Hydrogen s-wave reduced radial functions (n=1,2,3)"
set key top right
set grid
set xlabel "r"
set ylabel "f(r)"
plot \
  "wave_task_n1.dat" using 1:2 with lines lw 2 title "n=1 numeric", \
  "wave_task_n1.dat" using 1:3 with lines lw 2 dt 2 title "n=1 analytical", \
  "wave_task_n2.dat" using 1:2 with lines lw 2 title "n=2 numeric", \
  "wave_task_n2.dat" using 1:3 with lines lw 2 dt 2 title "n=2 analytical", \
  "wave_task_n3.dat" using 1:2 with lines lw 2 title "n=3 numeric", \
  "wave_task_n3.dat" using 1:3 with lines lw 2 dt 2 title "n=3 analytical"
