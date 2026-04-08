CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
TXTFILES = sho.data damped_oscillator.data lotka_volterra.data orbit_circular.data orbit_newtonian_elliptic.data orbit_relativistic_precession.data three_body_figure8.data

all : plot



Out.txt ${TXTFILES}: main
	./main > Out.txt


plot: plots.svg orbit_plot figure8_plot

plots.svg: plot.gnu ${TXTFILES}
	gnuplot plot.gnu

orbit_plot: orbit_plots.svg

orbit_plots.svg: plot_orbit.gnu ${TXTFILES}
	gnuplot plot_orbit.gnu

figure8_plot: figure8_plot.svg

figure8_plot.svg: plot_three_body.gnu ${TXTFILES}
	gnuplot plot_three_body.gnu

main: main.o ode_solver.o vector.o

main.o: main.cc ode_solver.h vector.h

ode_solver.o: ode_solver.cc ode_solver.h vector.h

vector.o: vector.cc vector.h




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data