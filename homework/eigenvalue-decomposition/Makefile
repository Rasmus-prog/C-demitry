CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

RMAX ?= 10
DR ?= 0.3
DATFILES = wavefunction_n1.dat wavefunction_n2.dat wavefunction_n3.dat


all : Out.txt dat

Out.txt: main
	./main --jacobi-test > $@


dat: $(DATFILES)

$(DATFILES): main
	./main -rmax $(RMAX) -dr $(DR) > hydrogen_run.txt

main: main.o vector.o matrix.o evd.o

main.o: vector.h matrix.h evd.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

evd.o: vector.h matrix.h evd.h evd.cc



.PHONY : clean dat
clean:
	$(RM) [Oo]ut* main *.o *.dat