CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm


all: spline_quadratic.svg spline_linear.svg spline_cubic.svg spline_cubic_check.svg data.dat


data.dat: main
	./main > $@


spline_quadratic.svg: data.dat spline_quadratic.gnu
	gnuplot spline_quadratic.gnu


spline_linear.svg: data.dat spline_linear.gnu
	gnuplot spline_linear.gnu


spline_cubic.svg: data.dat spline_cubic.gnu
	gnuplot spline_cubic.gnu


spline_cubic_check.svg: data.dat spline_cubic_check.gnu
	gnuplot spline_cubic_check.gnu


main: main.o vector.o qspline.o cspline.o

main.o: vector.h qspline.h cspline.h main.cc

cspline.o: cspline.h vector.h cspline.cc

qspline.o: qspline.h vector.h qspline.cc

vector.o: vector.h vector.cc




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.dat