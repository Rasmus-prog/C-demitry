# gnuplot script for visualizing QR decomposition timings
# usage: gnuplot -persist plot.gnu
set terminal svg background rgb "#ffffff"
set output "qr_decomposition_runtime.svg"
set title "QR decomposition runtime vs matrix size"
set xlabel "Matrix size N"
set ylabel "Elapsed time (s)"
set key left top
set grid

# cubic fit
f(x) = a * x**3
fit f(x) "out.times.data" using 1:2 via a

# plot data points and fitted curve
plot "out.times.data" using 1:2 with points pt 7 title "measured", \
     f(x) title sprintf("fit: %g*N^3", a)
