set terminal svg background rgb "white"
set output "scaling.svg"
set title "Monte Carlo Error Scaling"
set xlabel "Number of Samples"
set ylabel "Actual Error"
set logscale x
set logscale y
set grid

plot"scaling.data" using 1:4 with linespoints title "std::mt19937", \
    "scaling.data" using 1:5 with linespoints title "Quasi-random", \
    "scaling.data" using 1:6 with linespoints title "LCG", \
     1/sqrt(x) with lines title "1/√N"