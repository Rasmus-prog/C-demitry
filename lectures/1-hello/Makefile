CXXFLAFS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

Out.txt : main
		./main > Out.txt

main : main.o hello.o
		$(CXX) main.o hello.o -o main $(LDLIBS)

main.o : hello.h main.cc
		$(CXX) $(CXXFLAFS) -c main.cc -o main.o

hello.o : hello.h hello.cc
		$(CXX) $(CXXFLAFS) -c hello.cc -o hello.o


.PHONEY:clean
clean:
		$(RM) main main.o Out.txt hello.o

.PHONEY:test
test:
		@echo $(CXX)
		@echo CXX
		@echo $(RM) 