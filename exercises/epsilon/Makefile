CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

all : Out.txt
	@cat Out.txt


Out.txt: main
	./main > Out.txt

main.o : main.cc
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o
	$(CXX) -o main main.o $(LDLIBS)

.PHONY : clean
clean :
	$(RM) -f *.o main Out.txt