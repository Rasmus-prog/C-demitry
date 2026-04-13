CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
PYTHON ?= python3


all: plots Out.txt

plots : erf_acc.svg erf_out.svg error_quality.svg

error_quality.svg : error_quality.gnu error_quality.data
	gnuplot error_quality.gnu

erf_out.svg : erf_out.gnu erf_out.data
	gnuplot erf_out.gnu

erf_acc.svg : erf_acc.gnu erf_acc.data
	gnuplot erf_acc.gnu

erf_acc.data erf_out.data error_quality.data: main
	./main > /dev/null

Out.txt: main integration.py
	./main > $@
	printf '\n--- integration.py output ---\n' >> $@
	$(PYTHON) integration.py >> $@

main: main.o sfuns.o

main.o: main.cc

sfuns.o: sfuns.cc sfuns.h


.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o *.svg *.data *.txt