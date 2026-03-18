set terminal svg background rgb "#ffffff"
set output 'spline_cubic.svg'
set xlabel 'x'
set ylabel 'value'
set title 'Cubic spline for cos(x) and its anti-derivative on [0,9]'
set key bottom right
set grid
set label 1 'Table: x_i = 0, 0.5, ..., 9 ; y_i = cos(x_i)' at graph 0.24, 0.92
plot 'data.dat' index 0 using 1:4 with lines lw 3 lc rgb '#d62728' title 'Cubic interpolant',\
     'data.dat' index 0 using 1:12 with lines lw 3 lc rgb '#1f77b4' title 'Interpolant anti-derivative',\
     'data.dat' index 0 using 1:5 with lines lw 1 dt 2 lc rgb '#2ca02c' title 'cos(x)',\
     'data.dat' index 0 using 1:13 with lines lw 1 dt 2 lc rgb '#ff7f0e' title 'sin(x)',\
     'data.dat' index 1 using 1:2 with points pt 7 ps 1.0 lc rgb '#000000' title 'Tabulated knots'
