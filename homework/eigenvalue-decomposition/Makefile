CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
SHELL := /bin/bash
RMAX ?= 10
DR ?= 0.3
DRS ?= 0.60 0.50 0.40 0.30 0.25 0.20 0.15 0.12 0.10
RMAXS ?= 4 6 8 10 12 15 18 22 26 30
NS ?= $(shell seq 20 5 70)

OBJS = main.o matrix.o jacobi.o


all: full

full: Out.txt plots convergence-dr-serial.txt convergence-dr-parallel.txt convergence-rmax-parallel.txt scaling-plots


main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDLIBS) -o $@

main.o: main.cc jacobi.h matrix.h
matrix.o: matrix.cc matrix.h
jacobi.o: jacobi.cc jacobi.h matrix.h

Out.txt: main
	./main -rmax $(RMAX) -dr $(DR) > $@

hydrogen-data: Out.txt

plots: hydrogen-data hydrogen_plots.gnu
	gnuplot hydrogen_plots.gnu

convergence-dr-serial.txt: main Makefile
	@start=$$(date +%s%N); \
	>$@; \
	for DRV in $(DRS); do \
		./main -rmax $(RMAX) -dr $$DRV -ground >> $@; \
	done; \
	sort -g -k2,2 $@ -o $@; \
	end=$$(date +%s%N); \
	echo $$(( (end-start)/1000000 )) > convergence-dr-serial.wall-ms.txt

convergence-dr-parallel.txt: main Makefile
	@start=$$(date +%s%N); \
	$(RM) log.dr.*; \
	for DRV in $(DRS); do \
		./main -rmax $(RMAX) -dr $$DRV -ground > log.dr.$$DRV & \
	done; \
	wait; \
	>$@; \
	for DRV in $(DRS); do cat log.dr.$$DRV >> $@; done; \
	sort -g -k2,2 $@ -o $@; \
	$(RM) log.dr.*; \
	end=$$(date +%s%N); \
	echo $$(( (end-start)/1000000 )) > convergence-dr-parallel.wall-ms.txt

convergence-rmax-parallel.txt: main Makefile
	@$(RM) log.rmax.*; \
	for RMAXV in $(RMAXS); do \
		./main -rmax $$RMAXV -dr $(DR) -ground > log.rmax.$$RMAXV & \
	done; \
	wait; \
	>$@; \
	for RMAXV in $(RMAXS); do cat log.rmax.$$RMAXV >> $@; done; \
	sort -g -k1,1 $@ -o $@; \
	$(RM) log.rmax.*

scaling.times.txt: main Makefile
	@$(RM) log.scale.*; \
	for N in $(NS); do \
		./main -benchmark $$N > log.scale.$$N & \
	done; \
	wait; \
	>$@; \
	for N in $(NS); do awk '{print $$1, $$2}' log.scale.$$N >> $@; done; \
	sort -g -k1,1 $@ -o $@; \
	$(RM) log.scale.*

scaling.ratio.txt: scaling.times.txt
	@awk '{printf "%d %.12g %.12g\n", $$1, $$2, $$2/($$1*$$1*$$1)}' $< > $@

scaling: scaling.times.txt scaling.ratio.txt

scaling-plots: scaling scaling_plots.gnu
	gnuplot scaling_plots.gnu



.PHONY : all full make clean hydrogen-data plots scaling scaling-plots convergence-dr-serial.txt convergence-dr-parallel.txt convergence-rmax-parallel.txt
clean:
	$(RM) [Oo]ut* main *.o hydrogen_*.txt *.svg convergence-*.txt convergence-*.wall-ms.txt scaling.times.txt scaling.ratio.txt log.*