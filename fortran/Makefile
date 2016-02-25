#$preamble
# A simple hand-made makefile nsga2

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
FCFLAGS = -O3
# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
# LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = assign_nsga2
# nsga2_mc_assign


# "make" builds all
all: $(PROGRAMS)

#$part2
# prog2.o: prog2.incf

#$part3
# prog3: aux.o

#$part4
# Using Fortran MODULES:
assign_nsga2.o: qsort_c_module.o nsga2_module.o assign_module.o
assign_nsga2: qsort_c_module.o nsga2_module.o assign_module.o

# nsga2_mc_assign.o: qsort_c_module.o nsga2_module.o
# nsga2_mc_assign: qsort_c_module.o nsga2_module.o


# General rules
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules.
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
