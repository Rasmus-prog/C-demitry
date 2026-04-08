CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm


all: plots

plots : erf_acc.svg erf_out.svg

erf_out.svg : erf_out.gnu erf_out.data
	gnuplot erf_out.gnu

erf_acc.svg : erf_acc.gnu erf_acc.data
	gnuplot erf_acc.gnu

erf_acc.data erf_out.data: main
	./main > Out.txt


Out.txt: main
	./main > $@

main: main.o

main.o: main.cc


.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data