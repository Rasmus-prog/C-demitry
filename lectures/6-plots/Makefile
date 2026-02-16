CXXFLAGS += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

all: fig1.svg fig2.png

test.txt: erf.data Makefile
	>$@
	cat erf.data | while read a b c; do echo $$a $$c >> $@; done

fig2.png: erf.data Makefile
	echo '\
	set terminal png; \
	set output "$@"; \
	set key right bottom; \
	set tics out; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot "$<" index 1 using 1:3 with lines \
	' | tee log.fig2.gnu | gnuplot 




fig1.svg: fig1.gnu erf.data
	gnuplot $<

erf.data: main Makefile
	./$< -xmin -2 -xmax 2 -dx 0.5 > $@ 2>log.erf
	echo "\n\n" >> $@
	./$< -xmin 0 -xmax 10 -dx 0.255 >> $@ 2>>log.erf


.PHONY: clean
clean:
	rm -f fig1.svg erf.data log.erf main fig2.png log.fig2.gnu test.txt