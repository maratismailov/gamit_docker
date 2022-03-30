# Makefile for Ghandlers on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:30:34 UTC 2021

all: cf_hndl_lib



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

CF_HNDL_LIB = cf_hndl_lib.a
SLIB = 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

cf_hndl_lib : $(CF_HNDL_LIB)

$(CF_HNDL_LIB) : \
	$(CF_HNDL_LIB)(check_cflag.o) \
	$(CF_HNDL_LIB)(check_csize.o) \
	$(CF_HNDL_LIB)(open_cf.o) \
	$(CF_HNDL_LIB)(rd_cf_actv.o) \
	$(CF_HNDL_LIB)(wr_cf_actv.o) \
	$(CF_HNDL_LIB)(create_cf.o) \
	$(CF_HNDL_LIB)(rw_gobs_util.o) \
	$(CF_HNDL_LIB)(mfile_util.o)
	ranlib cf_hndl_lib.a

$(CF_HNDL_LIB)(rd_cf_actv.o): ../includes/cfile_def.h ../includes/const_param.h
$(CF_HNDL_LIB)(wr_cf_actv.o): ../includes/cfile_def.h
$(CF_HNDL_LIB)(rw_gobs_util.o): ../includes/gobs_def.h ../includes/kalman_param.h
$(CF_HNDL_LIB)(mfile_util.o): ../includes/kalman_param.h ../includes/mfile_def.h

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
