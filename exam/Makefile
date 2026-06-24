CXXFLAGS += -std=c++23 -Wall -Wextra -O2
LDLIBS   += -lstdc++ -lm



all: run plots 

plots: comparison.svg sweep.svg

comparison.svg: comparison.gnu $(wildcard data/*.data)
	mkdir -p plots
	gnuplot comparison.gnu

sweep.svg: sweep.gnu $(wildcard data/sweep_*.data)
	mkdir -p plots
	gnuplot sweep.gnu

data/sqrt_x.data data/inv_sqrt_x.data data/log_sqrt_x.data data/circle.data \
data/sweep_sqrt_x.data data/sweep_inv_sqrt_x.data data/sweep_log_sqrt_x.data: main
	mkdir -p data
	./main



main: main.o sfuns.o quad_random.o

main.o:       main.cc sfuns.h quad_random.h
sfuns.o:      sfuns.cc sfuns.h
quad_random.o: quad_random.cc quad_random.h

.PHONY: clean run
run: main
	mkdir -p data && ./main

clean:
	$(RM) -r Out* main *.o data/ plots/