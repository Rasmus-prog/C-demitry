CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS += -lstdc++ -lm
SHELL := /bin/bash
RMAX ?= 10
DR ?= 0.3
DRS ?= 0.60 0.50 0.40 0.30 0.25 0.20 0.15 0.12 0.10
RMAXS ?= 4 6 8 10 12 15 18 22 26 30
NS ?= $(shell seq 20 5 70)

OBJS = main.o matrix.o jacobi.o

# Default target
all: hydrogen scaling-plots clean-logs

# Main executable
main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDLIBS) -o $@

main.o: main.cc jacobi.h matrix.h
matrix.o: matrix.cc matrix.h
jacobi.o: jacobi.cc jacobi.h matrix.h

# ============================================
# HYDROGEN ATOM, S-WAVE RADIAL EQUATION
# ============================================

hydrogen: Out.txt hydrogen-plots hydrogen-data

Out.txt: main
	@echo "Computing hydrogen atom solution..."
	@./main -full -rmax $(RMAX) -dr $(DR) > $@

hydrogen-plots: hydrogen-data hydrogen_plots.gnu
	@echo "Generating hydrogen plots..."
	@gnuplot hydrogen_plots.gnu

hydrogen-data: main
	@./main -rmax $(RMAX) -dr $(DR) > /dev/null

# ============================================
# CONVERGENCE STUDIES (Serial vs Parallel)
# ============================================

convergence: convergence-dr-serial.txt convergence-dr-parallel.txt \
             convergence-rmax-parallel.txt

convergence-dr-serial.txt: main Makefile
	@echo "Computing Δr convergence (SERIAL)..."
	@start=$$(date +%s%N); \
	>$@; \
	for DRV in $(DRS); do \
		./main -rmax $(RMAX) -dr $$DRV -ground >> $@; \
	done; \
	sort -g -k2,2 $@ -o $@; \
	end=$$(date +%s%N); \
	echo $$(( (end-start)/1000000 )) > convergence-dr-serial.wall-ms.txt

convergence-dr-parallel.txt: main Makefile
	@echo "Computing Δr convergence (PARALLEL)..."
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
	@echo "Computing rmax convergence (PARALLEL)..."
	@$(RM) log.rmax.*; \
	for RMAXV in $(RMAXS); do \
		./main -rmax $$RMAXV -dr $(DR) -ground > log.rmax.$$RMAXV & \
	done; \
	wait; \
	>$@; \
	for RMAXV in $(RMAXS); do cat log.rmax.$$RMAXV >> $@; done; \
	sort -g -k1,1 $@ -o $@; \
	$(RM) log.rmax.*

# ============================================
# SCALING & OPTIMIZATION
# ============================================

scaling: scaling-times scaling-ratio

scaling-times: scaling.times.txt

scaling-ratio: scaling.ratio.txt

scaling.times.txt: main Makefile
	@echo "Measuring O(n³) scaling (PARALLEL)..."
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

scaling-plots: scaling scaling_plots.gnu
	@echo "Generating performance plots..."
	@gnuplot scaling_plots.gnu

# ============================================
# VERIFICATION & TESTING
# ============================================

verify:
	@echo "Running verification tests..."
	@./main -verify 5
	@./main -verify 10

benchmark-opt:
	@echo "Running optimization comparison..."
	@./main -compare-opt

# ============================================
# HOUSEKEEPING
# ============================================

clean-logs:
	@$(RM) log.*

clean: clean-logs
	$(RM) [Oo]ut* main *.o hydrogen_*.txt *.svg convergence-*.txt \
	      convergence-*.wall-ms.txt scaling.times.txt scaling.ratio.txt

# Phony targets
.PHONY: all hydrogen convergence scaling scaling-times scaling-ratio scaling-plots \
        verify benchmark-opt clean clean-logs hydrogen-plots hydrogen-data