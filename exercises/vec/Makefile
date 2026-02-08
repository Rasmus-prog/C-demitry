CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

Out.txt: main
	./main > Out.txt

main.o : main.cc
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

vec.o : vec.cc
	$(CXX) $(CXXFLAGS) -o vec.o -c vec.cc

main : main.o vec.o
	$(CXX) -o main main.o vec.o $(LDLIBS)

.PHONY : clean
clean :
	$(RM) -f *.o main Out.txt