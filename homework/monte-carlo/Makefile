CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
PYTHON ?= python3


all: Out.txt plots

plots: unit_circle.svg


unit_circle.svg: unit_circle.gnu unit_circle.data
	gnuplot $<


Out.txt unit_circle.data: main 
	./main > $@

main: main.o sfuns.o

sfuns.o: sfuns.cc sfuns.h

main.o: main.cc 

sfuns.cc: sfuns.h



.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data *.txt