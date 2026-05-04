set terminal svg background rgb 'white'
set output 'plot.svg'
set datafile commentschars '#'

set title 'ANN interpolation of tabulated function'
set xlabel 'x'
set ylabel 'y'
set grid
set key left top

plot \
	'Out.txt' using 1:2 with lines lw 2 lc rgb '#1f77b4' title 'target function', \
	'Out.txt' using 1:3 with lines lw 2 dt 2 lc rgb '#d62728' title 'network response'
