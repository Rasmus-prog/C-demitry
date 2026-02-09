CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm

all: Out.txt Out-cmdline.txt

Out-cmdline.txt : main-cmdline minusn.txt data.txt
	./main-cmdline -dr 0.05 -rmax 20 -input data.txt $$(cat minusn.txt) > $@

Out.txt Err.txt: script.py main data.txt Makefile
	python3 $< 2> Err.txt
	./main < data.txt > Out.txt 2>> Err.txt

data.txt: Makefile
	echo 1 2 > $@
	echo 2 1 >> $@
	echo 3 3 >> $@
	echo 4 5 >> $@
	echo 5 6 >> $@

N=100
minusn.txt : Makefile
	>$@
	for i in `seq $(N)`; do echo "-n " $$i >> $@; done
.Phony : clean
clean :
		rm -f main Err.txt Out.txt out.txt *.o main-cmdline