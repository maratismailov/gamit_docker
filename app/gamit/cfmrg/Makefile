# Makefile for cfmrg on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:16:53 UTC 2021

all: cfmrg



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

CFMRG_LIB = cfmrg_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

cfmrg : cfmrg.f $(CFMRG_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) cfmrg.f $(CFMRG_LIB) $(SLIB) -o cfmrg
	rm -f cfmrg.o

$(CFMRG_LIB) : \
	$(CFMRG_LIB)(aprior.o) \
	$(CFMRG_LIB)(check_limits.o) \
	$(CFMRG_LIB)(cversn.o) \
	$(CFMRG_LIB)(fills1.o) \
	$(CFMRG_LIB)(fills2.o) \
	$(CFMRG_LIB)(getcfs.o) \
	$(CFMRG_LIB)(queryb.o) \
	$(CFMRG_LIB)(queryp.o)
	ranlib cfmrg_lib.a

$(CFMRG_LIB)(aprior.o): ../includes/dimpar.h
$(CFMRG_LIB)(check_limits.o): ../includes/dimpar.h
$(CFMRG_LIB)(fills1.o): ../includes/dimpar.h
$(CFMRG_LIB)(fills2.o): ../includes/dimpar.h
$(CFMRG_LIB)(getcfs.o): ../includes/dimpar.h
$(CFMRG_LIB)(queryb.o): ../includes/dimpar.h
$(CFMRG_LIB)(queryp.o): ../includes/dimpar.h

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
