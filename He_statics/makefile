# Makefile for program '4hedft' 
# 3D static (imaginary time evolution) 4-He Density Functional Theory
#

COMP = gfortran
CFLAGS = -c -O3 -fopenmp -fopt-info-vec-optimized -falign-commons -march=native -fcheck=all -fmax-errors=5 -ffree-line-length-512 -J./modules
LD_FLAGS = -lfftw3 -lm -fopenmp -lfftw3_threads

# Name of the program
PROGNAME = 4hedft

# Fortran objects
OBJS = init_deriv_parallel.o V_impur.o modules.o FT_V_spline.o BCN4HeDFT.o \
       derden.o	spline.o dimen.o energy.o evolo.o evolox.o \
       fforma.o fft.o initcg.o ironing.o mates.o \
       poten.o printout.o r_cm.o readen.o respar.o \
       term_alfa.o timer.o titols.o vareps.o varmu.o \
       tstgrid.o s13adf.o newder.o redef.o steprkr.o \
       steppcr.o steprkrl.o steppcrl.o diag.o instates.o \
       rhoasin0.o instates_external.o

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
