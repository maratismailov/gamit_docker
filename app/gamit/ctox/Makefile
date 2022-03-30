# Makefile for ctox on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:18:23 UTC 2021

all: ctox



# ASSIGMENTS

# ASSIGMENTS
### LF95 flags ###
#FTN = lf95
#FFLAGS = -O
#CC = gcc
#CFLAGS = -g
### IFC 7.0 flags ###
#FTN = ifc
#FFLAGS = -O -Vaxlib -save -zero -u -72 -w95 -w90 -cm
#CC = icc
#CFLAGS = -O
### IFORT 8.0 flags ###
#FTN = ifort
#FFLAGS = -O3 -Vaxlib -save -zero -u -72 -w95 -w90 -cm -assume byterecl -static -xN
#CC = icc
#CFLAGS = -O
### IFORT 9.0 flags ###
#FTN = ifort
#FFLAGS = -O3 -Vaxlib -save -zero -u -72 -w95 -w90 -cm -assume byterecl -static -Xt
#CC = icc
#CFLAGS = -O
### G77 flags ###
#FTN = g77
#FFLAGS = -O3 -Wuninitialized -Wunused -Wimplicit -fno-f2c -ffast-math -fno-automatic -fno-backslash -Wno-globals -fno-globals
#CC = gcc
#CFLAGS = -g
### GCC V4 gfortran flags ###
FTN = gfortran
FFLAGS = -O3 -Wuninitialized -fno-f2c -ffast-math -fno-automatic -fno-backslash -m64 -mcmodel=large
CC = gcc
#CFLAGS = -g
#### For 64-bit machines:
CFLAGS = -g -m64 -DADD64BIT

CTOX_LIB = ctox_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

ctox : ctox.f $(CTOX_LIB) $(SLIB) ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) ctox.f $(CTOX_LIB) $(SLIB) -o ctox
	rm -f ctox.o

$(CTOX_LIB) : \
	$(CTOX_LIB)(cdtred.o) \
	$(CTOX_LIB)(chdred.o) \
	$(CTOX_LIB)(cversn.o) \
	$(CTOX_LIB)(raddeg.o) \
	$(CTOX_LIB)(xdtrit.o) \
	$(CTOX_LIB)(xhdrit.o)
	ranlib ctox_lib.a

$(CTOX_LIB)(cdtred.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(CTOX_LIB)(chdred.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(CTOX_LIB)(xdtrit.o): ../includes/dimpar.h ../includes/model.h
$(CTOX_LIB)(xhdrit.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h

.c.a:
	$(CC) -c $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f90.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

.SUFFIXES: .f90
# End of Makefile
