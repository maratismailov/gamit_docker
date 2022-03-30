# Makefile for fixdrv on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:19:29 UTC 2021

all: fixdrv read_sestbl



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

FIXDRV_LIB = fixdrv_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

fixdrv : fixdrv.f $(FIXDRV_LIB) $(SLIB) ../includes/dimpar.h ../includes/global.h
	$(FTN) $(FFLAGS) fixdrv.f $(FIXDRV_LIB) $(SLIB) -o fixdrv
	rm -f fixdrv.o

read_sestbl : read_sestbl.f  $(SLIB) 
	$(FTN) $(FFLAGS) read_sestbl.f  $(SLIB) -o read_sestbl
	rm -f read_sestbl.o

$(FIXDRV_LIB) : \
	$(FIXDRV_LIB)(acmake.o) \
	$(FIXDRV_LIB)(armake.o) \
	$(FIXDRV_LIB)(bmake.o) \
	$(FIXDRV_LIB)(chdred.o) \
	$(FIXDRV_LIB)(chrnum.o) \
	$(FIXDRV_LIB)(clkprm.o) \
	$(FIXDRV_LIB)(cmmake.o) \
	$(FIXDRV_LIB)(dcheck.o) \
	$(FIXDRV_LIB)(dbmake.o) \
	$(FIXDRV_LIB)(fversn.o) \
	$(FIXDRV_LIB)(imakef.o) \
	$(FIXDRV_LIB)(lfntop.o) \
	$(FIXDRV_LIB)(mdmake.o) \
	$(FIXDRV_LIB)(newchr.o) \
	$(FIXDRV_LIB)(nextt1.o) \
	$(FIXDRV_LIB)(ngmake.o) \
	$(FIXDRV_LIB)(nxepoc.o) \
	$(FIXDRV_LIB)(nxsat.o) \
	$(FIXDRV_LIB)(rdsitt.o) \
	$(FIXDRV_LIB)(seschk.o) \
	$(FIXDRV_LIB)(somake.o) \
	$(FIXDRV_LIB)(tgsats.o) \
	$(FIXDRV_LIB)(timeup.o) \
	$(FIXDRV_LIB)(upchr.o) \
	$(FIXDRV_LIB)(xcheck.o) \
	$(FIXDRV_LIB)(xcsats.o)
	ranlib fixdrv_lib.a

$(FIXDRV_LIB)(acmake.o): ../includes/dimpar.h
$(FIXDRV_LIB)(armake.o): ../includes/dimpar.h
$(FIXDRV_LIB)(bmake.o): ../includes/dimpar.h ../includes/global.h
$(FIXDRV_LIB)(chdred.o): ../includes/dimpar.h
$(FIXDRV_LIB)(cmmake.o): ../includes/dimpar.h
$(FIXDRV_LIB)(dcheck.o): ../includes/dimpar.h
$(FIXDRV_LIB)(dbmake.o): ../includes/dimpar.h
$(FIXDRV_LIB)(imakef.o): ../includes/dimpar.h
$(FIXDRV_LIB)(ngmake.o): ../includes/dimpar.h
$(FIXDRV_LIB)(nxepoc.o): ../includes/dimpar.h
$(FIXDRV_LIB)(nxsat.o): ../includes/dimpar.h
$(FIXDRV_LIB)(seschk.o): ../includes/dimpar.h
$(FIXDRV_LIB)(somake.o): ../../gamit/includes/dimpar.h
$(FIXDRV_LIB)(tgsats.o): ../includes/dimpar.h
$(FIXDRV_LIB)(xcheck.o): ../includes/dimpar.h
$(FIXDRV_LIB)(xcsats.o): ../includes/dimpar.h

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
