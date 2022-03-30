# Makefile for svpos on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:39:38 UTC 2021

all: svpos svdiff svsp3 svpdiff svtrack svsnr fitsnr spcsnr simsnr fixms modear



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

MODEAR_LIB = modear_lib.a
SLIB = ../Khandlers/KalObs_hndl_lib.a ../Ghandlers/cf_hndl_lib.a ../track/track_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

svpos : svpos.f  $(SLIB) ../../libraries/includes/const_freq.h ../includes/const_param.h svpos.h
	$(FTN) $(FFLAGS) svpos.f  $(SLIB) -o svpos
	rm -f svpos.o

svdiff : svdiff.f  $(SLIB) ../includes/const_param.h svdiff.h
	$(FTN) $(FFLAGS) svdiff.f  $(SLIB) -o svdiff
	rm -f svdiff.o

svsp3 : svsp3.f  $(SLIB) ../includes/const_param.h svsp3.h
	$(FTN) $(FFLAGS) svsp3.f  $(SLIB) -o svsp3
	rm -f svsp3.o

svpdiff : svpdiff.f  $(SLIB) 
	$(FTN) $(FFLAGS) svpdiff.f  $(SLIB) -o svpdiff
	rm -f svpdiff.o

svtrack : svtrack.f  $(SLIB) ../includes/const_param.h svtrack.h
	$(FTN) $(FFLAGS) svtrack.f  $(SLIB) -o svtrack
	rm -f svtrack.o

svsnr : svsnr.f  $(SLIB) ../includes/const_param.h svsnr.h
	$(FTN) $(FFLAGS) svsnr.f  $(SLIB) -o svsnr
	rm -f svsnr.o

fitsnr : fitsnr.f  $(SLIB) ../../libraries/includes/const_freq.h ../includes/const_param.h fitsnr.h
	$(FTN) $(FFLAGS) fitsnr.f  $(SLIB) -o fitsnr
	rm -f fitsnr.o

spcsnr : spcsnr.f  $(SLIB) ../../libraries/includes/const_freq.h ../../libraries/includes/freq_def.h ../includes/cfile_def.h ../includes/const_param.h spcsnr.h
	$(FTN) $(FFLAGS) spcsnr.f  $(SLIB) -o spcsnr
	rm -f spcsnr.o

simsnr : simsnr.f  $(SLIB) ../includes/const_param.h simsnr.h
	$(FTN) $(FFLAGS) simsnr.f  $(SLIB) -o simsnr
	rm -f simsnr.o

fixms : fixms.f  $(SLIB) ../../libraries/includes/const_freq.h ../includes/const_param.h fixms.h
	$(FTN) $(FFLAGS) fixms.f  $(SLIB) -o fixms
	rm -f fixms.o

modear : modear.f $(MODEAR_LIB) $(SLIB) ../includes/cfile_def.h ../includes/kalman_param.h modear.h
	$(FTN) $(FFLAGS) modear.f $(MODEAR_LIB) $(SLIB) -o modear
	rm -f modear.o

$(MODEAR_LIB) : \
	$(MODEAR_LIB)(update_cfi.o) \
	$(MODEAR_LIB)(read_nav.o) \
	$(MODEAR_LIB)(preread_cf.o) \
	$(MODEAR_LIB)(comp_clock.o) \
	$(MODEAR_LIB)(earth_to_inert.o) \
	$(MODEAR_LIB)(phase_mod.o)
	ranlib modear_lib.a

$(MODEAR_LIB)(update_cfi.o): ../../libraries/includes/const_freq.h ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h modear.h
$(MODEAR_LIB)(read_nav.o): ../includes/const_param.h modear.h
$(MODEAR_LIB)(preread_cf.o): ../includes/cfile_def.h ../includes/kalman_param.h modear.h
$(MODEAR_LIB)(comp_clock.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h modear.h
$(MODEAR_LIB)(earth_to_inert.o): ../includes/cfile_def.h ../includes/const_param.h
$(MODEAR_LIB)(phase_mod.o): ../includes/const_param.h

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
