set terminal svg background rgb "#ffffff"
set output 'fit.svg'
set xlabel 'Time t (days)'
set ylabel 'Activity (relative units)'
set title 'Radioactive decay of ThX (Rutherford and Soddy, 1902)'
set key top right
set grid
set label 1 'Uncertainty combinations: (c0+-dc0, c1+-dc1)' at graph 0.43, 0.88
plot 'plot.data' using 1:2:3 with yerrorbars title 'Data' pt 7 ps 1.3 lc rgb '#1f77b4',\
     'plot.curves.data' using 1:2 with lines lw 3 lc rgb '#d62728' title 'Best fit',\
     'plot.curves.data' using 1:3 with lines lw 1 dt 2 lc rgb '#2ca02c' title '(+, +)',\
     'plot.curves.data' using 1:4 with lines lw 1 dt 2 lc rgb '#ff7f0e' title '(+, -)',\
     'plot.curves.data' using 1:5 with lines lw 1 dt 3 lc rgb '#9467bd' title '(-, +)',\
     'plot.curves.data' using 1:6 with lines lw 1 dt 3 lc rgb '#8c564b' title '(-, -)'
