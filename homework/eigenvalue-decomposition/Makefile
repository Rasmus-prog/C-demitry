CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

RMAX ?= 15
DR ?= 0.3
DATFILES = wavefunction_n1.dat wavefunction_n2.dat wavefunction_n3.dat
CONVFILES = convergence_dr.dat convergence_rmax.dat


all : Out.txt dat plot

Out.txt: main
	./main --jacobi-test > $@


dat: $(DATFILES)

convergence: $(CONVFILES)

plot: convergence.svg

convergence.svg: convergence.gnu $(CONVFILES)
	gnuplot convergence.gnu

$(DATFILES): main
	./main -rmax $(RMAX) -dr $(DR) > hydrogen_run.txt

$(CONVFILES): main
	./main --convergence --conv-fixed-rmax $(RMAX) --conv-fixed-dr $(DR) > convergence_run.txt

main: main.o vector.o matrix.o evd.o

main.o: vector.h matrix.h evd.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

evd.o: vector.h matrix.h evd.h evd.cc



.PHONY : clean dat convergence plot
clean:
	$(RM) [Oo]ut* main *.o *.dat *.svg *.txt