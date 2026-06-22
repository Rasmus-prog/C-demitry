CXXFLAGS += -std=c++23 -Wall -Wextra -O2
LDLIBS   += -lstdc++ -lm

Out.txt: main
	./main > $@

main: main.o sfuns.o quad_random.o

main.o:       main.cc sfuns.h quad_random.h
sfuns.o:      sfuns.cc sfuns.h
quad_random.o: quad_random.cc quad_random.h

.PHONY: clean
clean:
	$(RM) Out* main *.o
