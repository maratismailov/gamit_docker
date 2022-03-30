# Makefile for displace on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:40:42 UTC 2021

all: displace_lib displace simplexfd



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

DISPLACE_LIB = displace_lib.a
SLIB = ./displace_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

displace_lib : $(DISPLACE_LIB)

displace : displace.f  $(SLIB) param.displace.h
	$(FTN) $(FFLAGS) displace.f  $(SLIB) -o displace
	rm -f displace.o

simplexfd : simplexfd.c  $(SLIB) 
	$(CC) $(CFLAGS) simplexfd.c $(SLIB) -lm  -o simplexfd
	rm -f simplexfd.o

$(DISPLACE_LIB) : \
	$(DISPLACE_LIB)(conic.o) \
	$(DISPLACE_LIB)(linalg.o) \
	$(DISPLACE_LIB)(longlatxy.o) \
	$(DISPLACE_LIB)(nreq.o) \
	$(DISPLACE_LIB)(lambert2xy.o) \
	$(DISPLACE_LIB)(ran.o) \
	$(DISPLACE_LIB)(lower_cholesky_factor.o) \
	$(DISPLACE_LIB)(inverse_check.o) \
	$(DISPLACE_LIB)(invert_upper_triangle_matrix.o) \
	$(DISPLACE_LIB)(square_upper_triangle_matrix.o)
	ranlib displace_lib.a

$(DISPLACE_LIB)(lower_cholesky_factor.o): math_tools_for_matrix_algebra.h
$(DISPLACE_LIB)(inverse_check.o): math_tools_for_matrix_algebra.h
$(DISPLACE_LIB)(invert_upper_triangle_matrix.o): math_tools_for_matrix_algebra.h
$(DISPLACE_LIB)(square_upper_triangle_matrix.o): math_tools_for_matrix_algebra.h

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
