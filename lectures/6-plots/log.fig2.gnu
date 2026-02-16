\
set terminal png; \
set output "fig2.png"; \
set key right bottom; \
set tics out; \
set xlabel "x"; \
set ylabel "y"; \
plot "erf.data" index 1 using 1:3 with lines \

