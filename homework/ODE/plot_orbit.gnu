set terminal svg background rgb "#ffffff"
set output 'orbit_plots.svg'

set size square
set grid
set key top right box

set xlabel 'x = r cos(phi)'
set ylabel 'y = r sin(phi)'
set title 'Planetary Orbit Comparison from u(phi) = 1/r(phi)'
plot \
	'out_orbit_circular.txt' using ((1/$2)*cos($1)):((1/$2)*sin($1)) with lines lw 2.5 dt 1 lc rgb '#1f77b4' title 'Newtonian Circular (eps = 0, u''(0)=0)', \
	'out_orbit_newtonian_elliptic.txt' using ((1/$2)*cos($1)):((1/$2)*sin($1)) with lines lw 2.5 dt 2 lc rgb '#d62728' title 'Newtonian Elliptic (eps = 0, u''(0)=-0.5)', \
	'out_orbit_relativistic_precession.txt' using ((1/$2)*cos($1)):((1/$2)*sin($1)) with lines lw 2.5 dt 3 lc rgb '#2ca02c' title 'Relativistic Precession (eps = 0.01)'
