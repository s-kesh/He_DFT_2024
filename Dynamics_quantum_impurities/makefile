#
# Makefile for program '4hetddft-isotropic' 
# 3D dynamic (real-time evolution) 4-He Density Functional Theory for isotropic PES's
#

# Compiler and flags
COMP = gfortran
CFLAGS = -c -O3 -fopenmp -fopt-info-vec-optimized -falign-commons -march=native \
	 -fcheck=all -fmax-errors=5 -ffree-line-length-512 -J$(MOD_DIR)
LD_FLAGS = -lfftw3 -lm -fopenmp -lfftw3_threads

# Directories
MOD_DIR = ./modules

# Name of the program
PROGNAME = 4hetddft-isotropic_wf

# Fortran objects
OBJS = init_deriv_parallel.o modules.o V_impur.o FT_V_spline.o DFT4He3d.o derden.o \
       dimen.o 	energy.o fforma.o fft.o initcg.o \
       morse.o 	mates.o poten.o printoutc.o r_cm.o \
       redef.o 	readenc.o respar.o term_alfa.o timer.o \
       titols.o tstgrid.o s13adf.o newder.o steprk.o \
       steppc.o potenimp.o potenimpini.o rhoasin0.o

# Rules
$(PROGNAME): $(OBJS)
	$(COMP) -o $(PROGNAME) $(OBJS) $(LD_FLAGS)

%.o: %.f90
	$(COMP) $(CFLAGS) -o $@ $<

%.o: %.f
	$(COMP) $(CFLAGS) -o $@ $<

# Cleanup
clean:
	rm -f *.o *.bak *.lst $(MOD_DIR)/*.mod

distclean: clean
	rm -f $(PROGNAME)
