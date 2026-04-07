CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
TXTFILES = out_sho.txt out_damped_oscillator.txt out_lotka_volterra.txt out_orbit_circular.txt out_orbit_newtonian_elliptic.txt out_orbit_relativistic_precession.txt

all : plot



Out.txt ${TXTFILES}: main
	./main > Out.txt


plot: plots.svg orbit_plot

plots.svg: plot.gnu ${TXTFILES}
	gnuplot plot.gnu

orbit_plot: orbit_plots.svg

orbit_plots.svg: plot_orbit.gnu ${TXTFILES}
	gnuplot plot_orbit.gnu

main: main.o ode_solver.o vector.o

main.o: main.cc ode_solver.h vector.h

ode_solver.o: ode_solver.cc ode_solver.h vector.h

vector.o: vector.cc vector.h




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg