CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm


all: spline.svg spline_linear.svg data.dat


data.dat: main
	./main > $@


spline.svg: data.dat spline_quadratic.gnu
	gnuplot spline_quadratic.gnu


spline_linear.svg: data.dat spline_linear.gnu
	gnuplot spline_linear.gnu


main: main.o vector.o qspline.o

main.o: vector.h qspline.h main.cc

qspline.o: qspline.h vector.h qspline.cc

vector.o: vector.h vector.cc




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.dat