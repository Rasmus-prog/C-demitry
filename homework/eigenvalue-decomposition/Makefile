CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm




Out.txt: main
	./main > $@

main: main.o vector.o matrix.o evd.o

main.o: vector.h matrix.h evd.h main.cc

vector.o: vector.h vector.cc

matrix.o: vector.h matrix.h matrix.cc

evd.o: vector.h matrix.h evd.h evd.cc



.PHONY : clean
clean:
	$(RM) [Oo]ut* main *.o