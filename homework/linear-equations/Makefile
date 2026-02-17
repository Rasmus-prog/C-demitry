CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm
MAXERR = 1



all: Out.txt qr_decomposition_runtime.svg
Out.txt:main ; ./$< > $@
main:main.o matrix.o
main.o:matrix.h
matrix.o:matrix.h
clean:
	$(RM) [Oo]ut* main *.o qr_decomposition_runtime.svg fit.log

# timing data for QR decomposition on NxN matrices
out.times.data: main Makefile
	> $@
	for N in $(shell seq 100 10 200); do \
		/usr/bin/time --format "$$N %e" --output $@ --append \
			./main -size:$$N 1>/dev/null 2>/dev/null; \
	done


qr_decomposition_runtime.svg: plot.gnu out.times.data
	gnuplot $<

