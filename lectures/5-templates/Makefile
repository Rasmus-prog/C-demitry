CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm


all : Out.txt
	cat $<

Out.txt :main
	./$< > $@


.PHONY : clean
clean :
	rm -f main Out.txt

