CXXFLAGS += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

all: fig1.svg fig2.svg fig3.svg

fig1.svg: fig1.gnu erf.data erf-tab.data
	gnuplot $<

fig2.svg: fig2.gnu gamma.data gamma-tab.data
	gnuplot $<

fig3.svg: fig3.gnu lngamma.data lngamma-tab.data
	gnuplot $<


main: main.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)


erf.data: main Makefile
	./main > $@

gamma.data: main Makefile
	./main -xmin 0 -xmax 3 > $@

lngamma.data: main Makefile
	./main -xmin 0 -xmax 3 -log > $@


.PHONY: clean
clean:
	rm -f main erf.data gamma.data lngamma.data fig1.svg fig2.svg fig3.svg