CXXFLAFS = -Wall -Wextra -Werror -std=c++23 -o
LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

Out.txt : main
		./main > Out.txt

main : main.o hello.o sfuns.o
# 		$(CXX) main.o hello.o sfuns.o $(LDLIBS) -o main
		$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o main

main.o : hello.h sfuns.h main.cc
		$(CXX) $(CXXFLAFS) -c main.cc -o main.o

hello.o : hello.h hello.cc
		$(CXX) $(CXXFLAFS) -c hello.cc -o hello.o

sfuns.o : sfuns.h sfuns.cc
		$(CXX) $(CXXFLAFS) -c sfuns.cc -o sfuns.o


.PHONEY:clean
clean:
		$(RM) main main.o Out.txt *.o

.PHONEY:test
test:
		@echo $(CXX)
		@echo CXX
		@echo $(RM) 