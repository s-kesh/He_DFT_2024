# Makefile for 4hetddft-isotropic
# 3D DFT real-time evolution for 4-He with isotropic PES

# Compiler and flags
FC       := ifx
CFLAGS   := -c -O3 -xAVX -align array64byte -qopenmp -parallel \
            -qopt-matmul -unroll0 -module modules
LDFLAGS  := -qopenmp -qmkl=parallel -qopt-matmul
INCLUDES := -I$(MKLROOT)/include/fftw

# Directories
OBJDIR   := objects
MODDIR   := modules
BINDIR   := bin

# Program name
PROGNAME := 4hetddft-isotropic

# List of Fortran source files (no path)
SRC := \
	modules.f90	init_deriv_parallel.f90	spline.f90   V_impur.f90	DFT4He3d.f90	derden.f90 \
	dimen.f90	energy.f90	 fforma.f90	fft.f90		initcg.f90 \
	morse.f90	mates.f90	 poten.f90	printoutc.f90	r_cm.f90 \
	redef.f90	readenc.f90	 respar.f90	term_alfa.f90	timer.f90 \
	titols.f90	tstgrid.f90	 s13adf.f90	newder.f90 	steprk.f90 \
	steppc.f90	potenimp.f90 potenimpini.f90

# Object files in objects/
OBJ := $(SRC:%.f90=$(OBJDIR)/%.o)

# Default target
all: $(MODDIR) $(OBJDIR) $(PROGNAME)

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

# Clean targets
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod *.lst *.bak $(BINDIR)/*

distclean: clean
	rm -f $(PROGNAME)

.PHONY: all clean distclean
