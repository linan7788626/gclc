FITSLIB = $(HOME)/Projects/CMBtools/cfitsio/lib
FITSINC = $(HOME)/Projects/CMBtools/cfitsio/include
HPXDIR  = $(HOME)/Projects/CMBtools/healpix/Healpix_3.11
HPXINC  = $(HPXDIR)/include

CFLAGS = -O3 -heap-arrays -vec-report0 -mkl
FFLAGS = -L$(HPXDIR)/lib -lhealpix -L$(FITSLIB) -lcfitsio -O3 -mkl

INCLUDE = -I$(HPXINC) -I../utils

FC90 = mpif90

OBJ  = cmb.o gclc_sims.o

default: $(OBJ)

%.o: %.f90
	$(FC90) $(CFLAGS) $(INCLUDE) -c -o $*.o $*.f90

%.o: %.F90
	$(FC90) $(CFLAGS) $(INCLUDE) -c -o $*.o $*.F90

clean:
	rm -f *.mod *.o
