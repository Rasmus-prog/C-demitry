CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm

all : Out.txt

main: main.o sfuns.o

main.o: main.cc sfuns.h

sfuns.o: sfuns.cc sfuns.h

Out.txt: main
	./main > $@

.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o
