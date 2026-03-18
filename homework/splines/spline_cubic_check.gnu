set terminal svg background rgb "#ffffff"
set output 'spline_cubic_check.svg'
set xlabel 'x'
set ylabel 'value'
set title 'Cubic spline check: our cubic vs gnuplot built-in csplines'
set key bottom right
set grid
set label 1 'Black: our cubic spline (from data block 0), Blue dashed: gnuplot csplines (from knots)' at graph 0.02, 0.92
plot 'data.dat' index 0 using 1:4 with lines lw 3 lc rgb '#000000' title 'Our cubic spline',\
     'data.dat' index 1 using 1:2 smooth csplines with lines lw 2 dt 2 lc rgb '#1f77b4' title 'gnuplot smooth csplines',\
     'data.dat' index 1 using 1:2 with points pt 7 ps 0.9 lc rgb '#d62728' title 'Tabulated knots'
