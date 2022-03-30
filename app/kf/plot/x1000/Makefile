# Makefile for plot on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:40:59 UTC 2021

all: x1000



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

X1000 = x1000.a
SLIB = ../../../libraries/comlib/x1000/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

x1000 : $(X1000)

$(X1000) : \
	$(X1000)(check_window.o) \
	$(X1000)(cpux.o) \
	$(X1000)(cpuy.o) \
	$(X1000)(cwpy.o) \
	$(X1000)(cwpx.o) \
	$(X1000)(drawchar.o) \
	$(X1000)(drawmark.o) \
	$(X1000)(jiws.o) \
	$(X1000)(jmapm.o) \
	$(X1000)(jtexp.o) \
	$(X1000)(j2ply.o) \
	$(X1000)(mapem.o) \
	$(X1000)(nullterm.o) \
	$(X1000)(cufx_gftn.o) \
	$(X1000)(cufy_gftn.o) \
	$(X1000)(jfont_gftn.o) \
	$(X1000)(jbegn_gftn.o) \
	$(X1000)(j2drw_gftn.o) \
	$(X1000)(j2mov_gftn.o) \
	$(X1000)(jwind_gftn.o) \
	$(X1000)(jview_gftn.o) \
	$(X1000)(j2mrk_gftn.o) \
	$(X1000)(jcsiz_gftn.o) \
	$(X1000)(jclr_gftn.o) \
	$(X1000)(jcori_gftn.o) \
	$(X1000)(jjust_gftn.o) \
	$(X1000)(jmcur_gftn.o) \
	$(X1000)(jr2dr_gftn.o) \
	$(X1000)(jr2mv_gftn.o) \
	$(X1000)(jtexh_gftn.o) \
	$(X1000)(jwloc_gftn.o) \
	$(X1000)(jcolr_gftn.o) \
	$(X1000)(jlstl_gftn.o) \
	$(X1000)(jvtow_gftn.o) \
	$(X1000)(jend_gftn.o) \
	$(X1000)(gsclip_gftn.o)
	ranlib x1000.a


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
