CXXFLAGS += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

fig1.svg: fig1.gnu erf.data
	gnuplot $<

erf.data: main Makefile
	./$< -xmin -2 -xmax 2 dx 0.25 > $@ 2>log.erf


.PHONEY clean:
	rm -f fig1.svg erf.data log.erf main