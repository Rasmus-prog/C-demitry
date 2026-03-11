CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

all: Out.txt fit.svg

Out.txt: main
	./main > $@

fit.svg: Out.txt plot.gnu plot.data plot.curves.data
	gnuplot plot.gnu

main: main.o vector.o matrix.o qr.o lsfit.o

main.o: vector.h matrix.h qr.h lsfit.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

qr.o: vector.h matrix.h qr.h qr.cc

lsfit.o: vector.h matrix.h qr.h lsfit.h lsfit.cc

# # timing data for QR decomposition on NxN matrices
# out.times.data: main Makefile
# 	> $@
# 	for N in $(shell seq 100 20 400); do \
# 		/usr/bin/time --format "$$N %e" --output $@ --append \
# 			./main -size:$$N 1>/dev/null 2>/dev/null; \
# 	done


# qr_decomposition_runtime.svg: plot.gnu out.times.data
# 	gnuplot $<



.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.data *.svg