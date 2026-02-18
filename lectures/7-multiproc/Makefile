CXXFLAGS = -std=c++23 -Wall
LDLIBS += -lstdc++ -lm -lpthread
TIME = time --append --output $@ --format "$$i %e %U"


times.svg : Out.times Makefile
	@echo '\
		set terminal svg background rgb "white";\
		set output "$@";\
		set xlabel "Threads";\
		set ylabel "Time (seconds)";\
		a = 1;\
		b = 10;\
		f(n) = a+b/n;\
		fit f(x) "$<" using 1:2 via a,b;\
		plot \
			"$<" using 1:2 with linespoints title "Real Time", \
			 f(x) with lines title "Fitted Curve 1/n"\
		'| tee log.times | gnuplot


N = 1e9
Out.times : main Makefile
	>$@
	>log.main
	@for i in $$(seq 8); do  \
		$(TIME) ./main -terms $(N) -threads $$i  >> log.main; \
	done

all: out1 out2 out3 out4
out1 : main ; ./main -threads 1 -terms 1e8 > $@
out2 : main ; ./main -threads 1 -terms 2e8 > $@
out3 : main ; ./main -threads 1 -terms 3e8 > $@
out4 : main ; ./main -threads 1 -terms 4e8 > $@


test: main
	./main -threads 1 -terms 1e8 > out5 &
	./main -threads 1 -terms 2e8 > out6 &
	./main -threads 1 -terms 3e8 > out7 &
	./main -threads 1 -terms 4e8 > out8 &


.PHONEY: clean
clean:
	$(RM) -f log* [Oo]ut* main *.svg fit.log