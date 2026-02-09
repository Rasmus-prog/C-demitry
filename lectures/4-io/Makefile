CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm


Out.txt Err.txt : script.py main data.txt Makefile
	python3 $< 2> Err.txt
	./main < data.txt >> Out.txt 2>> Err.txt

data.txt : Makefile
	echo 1 2 > $@
	echo 2 1 >> $@
	echo 3 3 >> $@
	echo 4 5 >> $@
	echo 5 6 >> $@

.Phony : clean
clean :
		rm -f main Err.txt Out.txt out.txt *.o main-cmdline