CXX = g++ -std=c++23
CXXFLAGS = -Wall -Wextra #-Werror -O
LDLIBS += -lstdc++ -lm

all : Out.txt
	@cat Out.txt

Out.txt : main
# 		./main > Out.txt
		./$< > $@

main : main.o hello.o

main.o : main.cc hello.h vec3.h


.Phony : clean
clean :
		rm -f main Out.txt *.o