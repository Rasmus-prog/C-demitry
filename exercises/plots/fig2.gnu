set terminal svg background rgb "#ffffff"
set output "fig2.svg"
set xlabel "x"
set ylabel "Gamma(x)"
set key left top
set tics in
set grid xtics 
set grid ytics 
# set xrange [0:10]
plot \
    "gamma.data" using 1:3 with lines title "Gamma from sfuns", \
    "gamma-tab.data" using 1:2 with points pointtype 7 title "factorials"
