set terminal svg background rgb "#ffffff"
set output 'plots.svg'

set multiplot layout 3,1 title 'Embedded RK12 ODE Integrator Results'

set key left top
set grid
set xlabel 'x'
set ylabel 'u, u'''
plot \
    'out_sho.txt' using 1:2 with lines lw 2 title 'u(x) numeric', \
    'out_sho.txt' using 1:3 with lines lw 2 dt 2 title 'u''(x) numeric', \
    sin(x) with lines lw 1 dt 3 title 'sin(x) exact'

set xlabel 't'
set ylabel 'theta, omega'
plot \
    'out_damped_oscillator.txt' using 1:2 with lines lw 2 title 'theta(t)', \
    'out_damped_oscillator.txt' using 1:3 with lines lw 2 title 'omega(t)'

set xlabel 't'
set ylabel 'Population'
plot \
    'out_lotka_volterra.txt' using 1:2 with lines lw 2 title 'prey', \
    'out_lotka_volterra.txt' using 1:3 with lines lw 2 title 'predator'

unset multiplot
