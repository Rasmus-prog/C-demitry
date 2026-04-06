CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm




all : plot



Out.txt out_sho.txt out_damped_oscillator.txt out_lotka_volterra.txt: main
	./main > Out.txt


plot: plots.svg

plots.svg: plot.gnu out_sho.txt out_damped_oscillator.txt out_lotka_volterra.txt
	gnuplot plot.gnu

main: main.o vector.o

main.o: main.cc vector.h

vector.o: vector.cc vector.h




.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg