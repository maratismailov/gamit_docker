# Makefile for glorg on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:33:04 UTC 2021

all: glorg



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

GLORG_LIB = glorg_lib.a
SLIB = ../glout/glout_lib.a ../glfor/glfor_lib.a ../globk/globk_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glorg : glorg_main.f $(GLORG_LIB) $(SLIB) ../globk/globk_cmd_bd.f
	$(FTN) $(FFLAGS) glorg_main.f $(GLORG_LIB) $(SLIB) -o glorg
	rm -f glorg_main.o

$(GLORG_LIB) : \
	$(GLORG_LIB)(glorg.o) \
	$(GLORG_LIB)(glorg_util.o) \
	$(GLORG_LIB)(equate.o) \
	$(GLORG_LIB)(equate_parms.o) \
	$(GLORG_LIB)(decode_param.o) \
	$(GLORG_LIB)(force.o) \
	$(GLORG_LIB)(cov_xyz_neu.o) \
	$(GLORG_LIB)(apply_cfull.o) \
	$(GLORG_LIB)(est_plate.o) \
	$(GLORG_LIB)(merge_nonsec.o) \
	$(GLORG_LIB)(balance_cov.o)
	ranlib glorg_lib.a

$(GLORG_LIB)(glorg.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h ../includes/sd_common.h glorg_bd.f
$(GLORG_LIB)(glorg_util.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLORG_LIB)(equate.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h ../includes/orient.h
$(GLORG_LIB)(decode_param.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLORG_LIB)(force.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h ../includes/orient.h
$(GLORG_LIB)(cov_xyz_neu.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLORG_LIB)(est_plate.o): ../includes/const_param.h ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLORG_LIB)(merge_nonsec.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLORG_LIB)(balance_cov.o): ../includes/globk_common.h ../includes/kalman_param.h

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
