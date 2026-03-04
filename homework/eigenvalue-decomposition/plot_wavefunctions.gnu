set terminal svg background rgb "#ffffff"
set output "wavefunctions.svg"

set grid
set key outside right
set xlabel "r (Bohr)"
set ylabel "f(r)"
set title "Hydrogen s-wave wavefunctions"

set style line 1 lc rgb "#1f77b4" lw 2
set style line 2 lc rgb "#ff7f0e" lw 2
set style line 3 lc rgb "#2ca02c" lw 2

plot "wavefunctions.dat" using 1:2 with lines ls 1 title "n=1 num", \
     "wavefunctions.dat" using 1:3 with lines ls 1 dt 2 title "n=1 ex", \
     "wavefunctions.dat" using 1:4 with lines ls 2 title "n=2 num", \
     "wavefunctions.dat" using 1:5 with lines ls 2 dt 2 title "n=2 ex", \
     "wavefunctions.dat" using 1:6 with lines ls 3 title "n=3 num", \
     "wavefunctions.dat" using 1:7 with lines ls 3 dt 2 title "n=3 ex"

unset output
