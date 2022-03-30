# Makefile for blsum on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:36:38 UTC 2021

all: tslib blsum ensum tssum tsfit tscon tscomp xysum xyzsave enfit blavg bcsum makebase tsjumps enfit2apr eulvect



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

TSLIB = tslib.a
SLIB = ./tslib.a ../glfor/glfor_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

tslib : $(TSLIB)

blsum : blsum.f  $(SLIB) blsum.h
	$(FTN) $(FFLAGS) blsum.f  $(SLIB) -o blsum
	rm -f blsum.o

ensum : ensum.f  $(SLIB) ensum.h
	$(FTN) $(FFLAGS) ensum.f  $(SLIB) -o ensum
	rm -f ensum.o

tssum : tssum.f  $(SLIB) tssum.h
	$(FTN) $(FFLAGS) tssum.f  $(SLIB) -o tssum
	rm -f tssum.o

tsfit : tsfit.f  $(SLIB) ../includes/const_param.h tsfit.h
	$(FTN) $(FFLAGS) tsfit.f  $(SLIB) -o tsfit
	rm -f tsfit.o

tscon : tscon.f  $(SLIB) ../includes/const_param.h tscon.h tsfit.h tssum.h
	$(FTN) $(FFLAGS) tscon.f  $(SLIB) -o tscon
	rm -f tscon.o

tscomp : tscomp.f  $(SLIB) ../includes/const_param.h tscon.h tsfit.h
	$(FTN) $(FFLAGS) tscomp.f  $(SLIB) -o tscomp
	rm -f tscomp.o

xysum : xysum.f  $(SLIB) ensum.h
	$(FTN) $(FFLAGS) xysum.f  $(SLIB) -o xysum
	rm -f xysum.o

xyzsave : xyzsave.f  $(SLIB) ensum.h
	$(FTN) $(FFLAGS) xyzsave.f  $(SLIB) -o xyzsave
	rm -f xyzsave.o

enfit : enfit.f  $(SLIB) ../includes/const_param.h enfit.h
	$(FTN) $(FFLAGS) enfit.f  $(SLIB) -o enfit
	rm -f enfit.o

blavg : blavg.f  $(SLIB) blavg.h
	$(FTN) $(FFLAGS) blavg.f  $(SLIB) -o blavg
	rm -f blavg.o

bcsum : bcsum.f  $(SLIB) bcsum.h
	$(FTN) $(FFLAGS) bcsum.f  $(SLIB) -o bcsum
	rm -f bcsum.o

makebase : makebase.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) makebase.f  $(SLIB) -o makebase
	rm -f makebase.o

tsjumps : tsjumps.f  $(SLIB) tsjumps.h
	$(FTN) $(FFLAGS) tsjumps.f  $(SLIB) -o tsjumps
	rm -f tsjumps.o

enfit2apr : enfit2apr.f  $(SLIB) 
	$(FTN) $(FFLAGS) enfit2apr.f  $(SLIB) -o enfit2apr
	rm -f enfit2apr.o

eulvect : eulvect.f  $(SLIB) 
	$(FTN) $(FFLAGS) eulvect.f  $(SLIB) -o eulvect
	rm -f eulvect.o

$(TSLIB) : \
	$(TSLIB)(tsutil.o) \
	$(TSLIB)(tsdesc.o) \
	$(TSLIB)(tsread.o) \
	$(TSLIB)(tsutcon.o) \
	$(TSLIB)(get_unr_ant.o)
	ranlib tslib.a

$(TSLIB)(tsutil.o): ../includes/const_param.h tsfit.h
$(TSLIB)(tsread.o): ../includes/const_param.h tssum.h
$(TSLIB)(tsutcon.o): ../includes/const_param.h tscon.h tsfit.h tssum.h
$(TSLIB)(get_unr_ant.o): tscon.h tssum.h

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
