CXX = g++
CXXFLAGS = -std=c++17 -O2

all: erf_plot

erf_plot: main.cc
	$(CXX) $(CXXFLAGS) main.cc -o erf_plot

clean:
	rm -f erf_plot tab.dat curve.dat
