CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

RMAX ?= 15
DR ?= 0.3
CONVFILES = convergence_dr.dat convergence_rmax.dat
WAVEFILES = wave_n1.dat wave_n2.dat wave_n3.dat


all : Out.txt plot waveplot

Out.txt: main
	./main --jacobi-test > $@


convergence: $(CONVFILES)

wave: $(WAVEFILES)

plot: convergence.svg

waveplot: wavefunctions.svg

convergence.svg: convergence.gnu $(CONVFILES)
	gnuplot convergence.gnu

wavefunctions.svg: wavefunctions.gnu $(WAVEFILES)
	gnuplot wavefunctions.gnu

$(CONVFILES): main
	./main --convergence --conv-fixed-rmax $(RMAX) --conv-fixed-dr $(DR) > convergence_run.txt

$(WAVEFILES): main
	./main --wavefunctions -rmax $(RMAX) -dr $(DR) --wf-states 3 > wavefunctions_run.txt

main: main.o vector.o matrix.o evd.o

main.o: vector.h matrix.h evd.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

evd.o: vector.h matrix.h evd.h evd.cc



.PHONY : clean convergence plot wave waveplot
clean:
	$(RM) [Oo]ut* main *.o *.dat *.svg *.txt