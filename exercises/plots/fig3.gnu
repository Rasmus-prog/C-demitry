set terminal svg background rgb "#ffffff"
set output "fig3.svg"
set xlabel "x"
set ylabel "lnGamma(x)"
set key left top
set tics in
set grid xtics 
set grid ytics 
# set xrange [0:10]
plot \
    "lngamma.data" using 1:3 with lines title "lnGamma from sfuns", \
    "lngamma-tab.data" using 1:2 with points pointtype 7 title "ln(n!)"
