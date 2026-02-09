CXXFLAGS = -Wall -Wextra -Werror -std=c++23
LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

Out.txt: main-cmd main-input-stream main-file
	./main > Out.txt

.PHONEY : clean
clean :
	$(RM) -f *.o main main-input main-file Out.txt out.txt outputfile.txt my_output_file.txt main-cmd main-input-stream main-file 