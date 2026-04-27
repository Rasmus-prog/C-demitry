CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
PYTHON ?= python3


all: Out.txt plots

plots: unit_circle.svg scaling.svg


unit_circle.svg: unit_circle.gnu unit_circle.data
	gnuplot $<

scaling.svg: scaling.gnu scaling.data
	gnuplot $<


Out.txt: main
	./main > $@

unit_circle.data: main
	./main > /dev/null

main: main.o sfuns.o lcg.o

sfuns.o: sfuns.cc sfuns.h

lcg.o: lcg.cc lcg.h

main.o: main.cc 

sfuns.cc: sfuns.h



.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data *.txt