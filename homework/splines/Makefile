CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm


all: spline.svg data.dat


data.dat: main
	./main > $@


spline.svg: data.dat spline.gnu
	gnuplot spline.gnu


main: main.o vector.o

main.o: vector.h main.cc

vector.o: vector.h vector.cc




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.dat