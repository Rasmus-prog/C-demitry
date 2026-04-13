set terminal svg background rgb "white"
set output "unit_circle.svg"
set title "Monte Carlo Estimation of Unit Circle error"
set xlabel "Number of Samples"
set ylabel "Estimated error"
set logscale x
set grid


plot "unit_circle.data" using 1:3 with linespoints title "Estimated Error", \
     "unit_circle.data" using 1:4 with linespoints title "Actual Error", \
     1/sqrt(x) with lines title "1/√N"

