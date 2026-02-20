CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm



all: qr_decomposition_runtime.svg Out.txt
Out.txt: main
	./main > $@

main: main.o matrix.o

main.o: matrix.h main.cc

matrix.o: matrix.h matrix.cc

# timing data for QR decomposition on NxN matrices
out.times.data: main Makefile
	> $@
	for N in $(shell seq 0 100 500); do \
		/usr/bin/time --format "$$N %e" --output $@ --append \
			./main -size:$$N 1>/dev/null 2>/dev/null; \
	done


qr_decomposition_runtime.svg: plot.gnu out.times.data
	gnuplot $<



.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o qr_decomposition_runtime.svg fit.log