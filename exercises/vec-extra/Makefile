CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm

Out.txt: main
	./main > Out.txt

main : main.cc vec.h
	$(CXX) $(CXXFLAGS) -o main main.cc $(LDLIBS)

.PHONY : clean
clean :
	$(RM) -f main Out.txt