#
# Makefile para montar el programa DFT4He3d 
#  (Density Functional Theory 4He 3dimensional program)
#
# (Version de compilacion para Pentium IV)
#
COMP =  ifort
CFLAGS = -c -O3 -xAVX -align array64byte -qopenmp\
		 -parallel -qopt-matmul -unroll0 -module ./modules
#CFLAGS = -c -g -check all -static -module ./modules 
LD_FLAGS = -threads -I${MKLROOT}/include/fftw -mkl=parallel -qopt-matmul 

#LD_FLAGS = -threads -parallel -qopt-matmul
#   Libraries for the FFT
#LIB1=fftw3
##LIB2=svml
#LIB3=fftw3_threads
#LIB4=pthread

#   Name of the program
PROGNAME=potential_multi_impurity


#   Fortran objects
objs=V_impur.o potential_multi_impurity.o spline.o
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
#	$(COMP)	-o $(PROGNAME) $(objs) -l$(LIB1) -l$(LIB3) -l$(LIB4)  $(LD_FLAGS)
#	$(COMP)	-o $(PROGNAME) $(objs) -l$(LIB1) -l$(LIB2)  $(LD_FLAGS)
	$(COMP)	-o $(PROGNAME) $(objs)  $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
modules:	
	$(COMP) $(CFLAGS) $(mods)

clear:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
clean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
