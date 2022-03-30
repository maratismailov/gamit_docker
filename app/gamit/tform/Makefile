# Makefile for tform on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:26:05 UTC 2021

all: tform conv_xyz2geo convertc



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

TFORM_LIB = tform_lib.a
CONV_XYZ2GEO_LIB = conv_xyz2geo_lib.a
CONVERTC_LIB = convertc_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

tform : tform.f $(TFORM_LIB) $(SLIB) ../includes/tform.h
	$(FTN) $(FFLAGS) tform.f $(TFORM_LIB) $(SLIB) -o tform
	rm -f tform.o

conv_xyz2geo : conv_xyz2geo.f $(CONV_XYZ2GEO_LIB) $(SLIB) ../includes/tform.h
	$(FTN) $(FFLAGS) conv_xyz2geo.f $(CONV_XYZ2GEO_LIB) $(SLIB) -o conv_xyz2geo
	rm -f conv_xyz2geo.o

convertc : convertc.f $(CONVERTC_LIB) $(SLIB) ../../libraries/includes/const_param.h
	$(FTN) $(FFLAGS) convertc.f $(CONVERTC_LIB) $(SLIB) -o convertc
	rm -f convertc.o

$(TFORM_LIB) : \
	$(TFORM_LIB)(carsph.o) \
	$(TFORM_LIB)(estimt.o) \
	$(TFORM_LIB)(estref.o) \
	$(TFORM_LIB)(forsys.o) \
	$(TFORM_LIB)(getcar.o) \
	$(TFORM_LIB)(getcor.o) \
	$(TFORM_LIB)(getcyl.o) \
	$(TFORM_LIB)(getgeo.o) \
	$(TFORM_LIB)(getloc.o) \
	$(TFORM_LIB)(getsname.o) \
	$(TFORM_LIB)(getsph.o) \
	$(TFORM_LIB)(gettfr.o) \
	$(TFORM_LIB)(putcar.o) \
	$(TFORM_LIB)(putcyl.o) \
	$(TFORM_LIB)(putgeo.o) \
	$(TFORM_LIB)(putloc.o) \
	$(TFORM_LIB)(putsph.o) \
	$(TFORM_LIB)(syminv.o) \
	$(TFORM_LIB)(tf2ref.o) \
	$(TFORM_LIB)(tversn.o)
	ranlib tform_lib.a

$(CONV_XYZ2GEO_LIB) : \
	$(CONV_XYZ2GEO_LIB)(putgeo.o) \
	$(CONV_XYZ2GEO_LIB)(getsname.o) \
	$(CONV_XYZ2GEO_LIB)(tversn.o)
	ranlib conv_xyz2geo_lib.a

$(CONVERTC_LIB) : \
	$(CONVERTC_LIB)(convertc_help.o) \
	$(CONVERTC_LIB)(carsph.o) \
	$(CONVERTC_LIB)(get_filetype.o) \
	$(CONVERTC_LIB)(read_glistf.o) \
	$(CONVERTC_LIB)(read_lfile.o) \
	$(CONVERTC_LIB)(read_velf.o) \
	$(CONVERTC_LIB)(write_lfile.o) \
	$(CONVERTC_LIB)(write_velf.o) \
	$(CONVERTC_LIB)(write_kmlf.o)
	ranlib convertc_lib.a

$(TFORM_LIB)(estimt.o): ../includes/tform.h
$(TFORM_LIB)(estref.o): ../includes/tform.h
$(TFORM_LIB)(forsys.o): ../includes/tform.h
$(TFORM_LIB)(getcar.o): ../includes/tform.h
$(TFORM_LIB)(getcor.o): ../includes/tform.h
$(TFORM_LIB)(getcyl.o): ../includes/tform.h
$(TFORM_LIB)(getgeo.o): ../includes/tform.h
$(TFORM_LIB)(getloc.o): ../includes/tform.h
$(TFORM_LIB)(getsph.o): ../includes/tform.h
$(TFORM_LIB)(gettfr.o): ../includes/tform.h
$(TFORM_LIB)(putcar.o): ../includes/tform.h
$(TFORM_LIB)(putcyl.o): ../includes/tform.h
$(TFORM_LIB)(putgeo.o): ../includes/tform.h
$(TFORM_LIB)(putloc.o): ../includes/tform.h
$(TFORM_LIB)(putsph.o): ../includes/tform.h
$(CONV_XYZ2GEO_LIB)(putgeo.o): ../includes/tform.h
$(CONVERTC_LIB)(read_glistf.o): ../../libraries/includes/const_param.h
$(CONVERTC_LIB)(read_lfile.o): ../../libraries/includes/const_param.h
$(CONVERTC_LIB)(read_velf.o): ../../libraries/includes/const_param.h
$(CONVERTC_LIB)(write_velf.o): ../../libraries/includes/const_param.h
$(CONVERTC_LIB)(write_kmlf.o): ../../libraries/includes/const_param.h

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
