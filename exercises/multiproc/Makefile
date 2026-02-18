CXXFLAGS = -std=c++23 -Wall
LDLIBS += -lstdc++ -lm -lpthread
TIME = time --append --output $@ --format "$$i %e %U"


times.svg : times.gnu Out.times Makefile
	gnuplot $<

N = 1e9
Out.times : main Makefile
	>$@
	>log.main
	@for i in $$(seq 10); do  \
		$(TIME) ./main -terms $(N) -threads $$i  >> log.main; \
	done


.PHONEY: clean
clean:
	$(RM) -f log* [Oo]ut* main *.svg fit.log