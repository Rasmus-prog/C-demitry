CXXFLAGS += -std=c++23 -Wall -Wextra -O2
LDLIBS += -lstdc++ -lm


all : Out.txt plots

plots : plot.gnu
	gnuplot $<

main: main.o sfuns.o

main.o: main.cc sfuns.h

sfuns.o: sfuns.cc sfuns.h

Out.txt: main
	./main > $@


.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data *.txt