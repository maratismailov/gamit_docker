# Makefile for matrix on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:12:51 UTC 2021

all: kinv_lib benchd



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

KINV_LIB = kinv_lib.a
SLIB = ./kinv_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

kinv_lib : $(KINV_LIB)

benchd : benchd.f  $(SLIB) 
	$(FTN) $(FFLAGS) benchd.f  $(SLIB) -o benchd
	rm -f benchd.o

$(KINV_LIB) : \
	$(KINV_LIB)(dwmab.o) \
	$(KINV_LIB)(dwadd.o) \
	$(KINV_LIB)(dvdot.o) \
	$(KINV_LIB)(dwdot.o) \
	$(KINV_LIB)(dwpiv.o) \
	$(KINV_LIB)(dwsmy.o) \
	$(KINV_LIB)(dwswp.o) \
	$(KINV_LIB)(dwsum.o) \
	$(KINV_LIB)(dwmov.o) \
	$(KINV_LIB)(dvswp.o) \
	$(KINV_LIB)(dwvmv.o) \
	$(KINV_LIB)(dwsad.o) \
	$(KINV_LIB)(dwsub.o) \
	$(KINV_LIB)(dwint.o) \
	$(KINV_LIB)(wmov.o) \
	$(KINV_LIB)(gauss_elim.o) \
	$(KINV_LIB)(invert_vis.o) \
	$(KINV_LIB)(scale_matrix.o) \
	$(KINV_LIB)(squez_matrix.o) \
	$(KINV_LIB)(dwadd8.o) \
	$(KINV_LIB)(dwint8.o) \
	$(KINV_LIB)(dwmov8.o) \
	$(KINV_LIB)(dwsmy8.o)
	ranlib kinv_lib.a


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
