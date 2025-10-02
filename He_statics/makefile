#
# Makefile for program '4hedft'
# 3D static (imaginary time evolution) 4-He Density Functional Theory
#

# Directories
OBJDIR   := objects
MODDIR   := modules
BINDIR   := bin

# Compiler
FC       := gfortran
INCLUDES := -I$(MODDIR) -I/usr/include

# Build type (default = release)
BUILD ?= release


# Flags for different builds
ifeq ($(BUILD),release)
  CFLAGS := -c -O3 -fopenmp -fopt-info-vec-optimized \
            -falign-commons -march=native -ffast-math \
            -fcheck=all -fmax-errors=5 -ffree-line-length-512 \
            -fomit-frame-pointer -fno-math-errno \
            -funroll-loops -J$(MODDIR)
  LDFLAGS := -fopenmp -lfftw3 -lm -lfftw3_threads
endif

ifeq ($(BUILD),debug)
  CFLAGS := -c -O0 -g -fcheck=all -Wall -Wextra -fbacktrace \
            -ffree-line-length-512 -J$(MODDIR)
  LDFLAGS := -g -fcheck=all -fbacktrace -lfftw3 -lm -lfftw3_threads
endif

ifeq ($(BUILD),profile)
  CFLAGS := -c -O2 -pg -fopenmp -march=native -ffast-math \
            -ffree-line-length-512 -J$(MODDIR)
  LDFLAGS := -pg -fopenmp -lfftw3 -lm -lfftw3_threads
endif

# Name of the program
PROGNAME = 4hedft

# List of Fortran source files
SRC := \
	modules.f90 init_deriv_parallel.f90 V_impur.f90 FT_V_spline.f90 \
	BCN4HeDFT.f90 derden.f90 spline.f90	dimen.f90 energy.f90 \
	evolo.f90 evolox.f90 fforma.f90	fft.f90 initcg.f90 ironing.f90 mates.f90 \
    poten.f90 printout.f90 r_cm.f90 readen.f90 respar.f90 \
    term_alfa.f90 timer.f90 titols.f90 vareps.f90 varmu.f90 \
    tstgrid.f90 s13adf.f90 newder.f90 redef.f90 steprkr.f90 \
    steppcr.f90 steprkrl.f90 steppcrl.f90 diag.f90 instates.f90 \
    rhoasin0.f90 instates_external.f90


# Object files in objects/
OBJ := $(SRC:%.f90=$(OBJDIR)/%.o)

# Default target
all: $(MODDIR) $(OBJDIR) $(BINDIR) $(PROGNAME)

# Link
$(PROGNAME): $(OBJ)
	$(FC) -o $(BINDIR)/$@ $^ $(INCLUDES) $(LDFLAGS)

# Compile Fortran 90 source
$(OBJDIR)/%.o: %.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(CFLAGS) -o $@ $<

# Compile Fortran 77 source
$(OBJDIR)/%.o: %.f | $(OBJDIR) $(MODDIR)
	$(FC) $(CFLAGS) -o $@ $<

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Clean targets
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod *.lst *.bak $(BINDIR)/*

distclean: clean
	rm -f $(PROGNAME)

.PHONY: all clean distclean
