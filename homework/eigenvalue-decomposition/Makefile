CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

RMAX ?= 15
DR ?= 0.1
CONVFILES = convergence_dr.dat convergence_rmax.dat
WAVEFILES = wave_n1.dat wave_n2.dat wave_n3.dat
TIMEFILE = out.times.data


all : Out.txt plot

Out.txt: main
	./main --jacobi-test > $@


convergence: $(CONVFILES)

wave: $(WAVEFILES)

plot: convergence.svg wavefunctions.svg timing.svg


convergence.svg: convergence.gnu $(CONVFILES)
	gnuplot convergence.gnu

wavefunctions.svg: wavefunctions.gnu $(WAVEFILES)
	gnuplot wavefunctions.gnu

timing.svg: timing.gnu $(TIMEFILE)
	gnuplot timing.gnu

$(CONVFILES): main
	./main --convergence --conv-fixed-rmax $(RMAX) --conv-fixed-dr $(DR) > convergence_run.txt

$(WAVEFILES): main
	./main --wavefunctions -rmax $(RMAX) -dr $(DR) --wf-states 3 > wavefunctions_run.txt

timing: $(TIMEFILE)

$(TIMEFILE): main Makefile
	> $@
	for N in $(shell seq 40 20 260); do \
		/usr/bin/time --format "$$N %e" --output $@ --append \
			./main --timing -size $$N 1>/dev/null 2>/dev/null; \
	done

main: main.o vector.o matrix.o evd.o

main.o: vector.h matrix.h evd.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

evd.o: vector.h matrix.h evd.h evd.cc



.PHONY : clean convergence plot wave waveplot timing
clean:
	$(RM) [Oo]ut* main *.o *.dat *.svg *.txt *.log