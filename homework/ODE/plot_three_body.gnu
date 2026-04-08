set terminal svg background rgb "#ffffff"
set output 'figure8_plot.svg'

set size square
set grid
set key top right box

set xlabel 'x'
set ylabel 'y'
set title 'Newtonian Three-Body Figure-8 Orbit (Equal Masses)'
plot \
	'three_body_figure8.data' using 8:9 with lines lw 8 lc rgb '#1f77b4' title 'Body 1', \
	'three_body_figure8.data' using 10:11 with lines lw 4 lc rgb '#d62728' title 'Body 2', \
	'three_body_figure8.data' using 12:13 with lines lw 2 lc rgb '#2ca02c' title 'Body 3'
