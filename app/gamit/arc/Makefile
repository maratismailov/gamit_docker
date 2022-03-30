# Makefile for arc on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:16:07 UTC 2021

all: arc



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

ARC_LIB = arc_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a ../lib/libsofa.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

arc : arc.f $(ARC_LIB) $(SLIB) ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
	$(FTN) $(FFLAGS) arc.f $(ARC_LIB) $(SLIB) -o arc
	rm -f arc.o

$(ARC_LIB) : \
	$(ARC_LIB)(adam.o) \
	$(ARC_LIB)(arcmrg.o) \
	$(ARC_LIB)(aversn.o) \
	$(ARC_LIB)(bilin8.o) \
	$(ARC_LIB)(binoml.o) \
	$(ARC_LIB)(calcof.o) \
	$(ARC_LIB)(check_gmodels.o) \
	$(ARC_LIB)(doodson_angle.o) \
	$(ARC_LIB)(dtwopi.o) \
	$(ARC_LIB)(earthrad.o) \
	$(ARC_LIB)(earthradTUM.o) \
	$(ARC_LIB)(eclout.o) \
	$(ARC_LIB)(eflux.o) \
	$(ARC_LIB)(ertorb.o) \
	$(ARC_LIB)(ERPFBOXW.o) \
	$(ARC_LIB)(eval.o) \
	$(ARC_LIB)(filopn.o) \
	$(ARC_LIB)(get_sat_info.o) \
	$(ARC_LIB)(gridval.o) \
	$(ARC_LIB)(hofrad.o) \
	$(ARC_LIB)(igs92.o) \
	$(ARC_LIB)(egm96.o) \
	$(ARC_LIB)(egm08.o) \
	$(ARC_LIB)(iers_etides.o) \
	$(ARC_LIB)(init.o) \
	$(ARC_LIB)(keplr.o) \
	$(ARC_LIB)(legnd2.o) \
	$(ARC_LIB)(legndr.o) \
	$(ARC_LIB)(legndr_p.o) \
	$(ARC_LIB)(matmpy.o) \
	$(ARC_LIB)(mattrn.o) \
	$(ARC_LIB)(oscrat.o) \
	$(ARC_LIB)(PROPBOXW.o) \
	$(ARC_LIB)(rdsrpgrd.o) \
	$(ARC_LIB)(read_arc_batch.o) \
	$(ARC_LIB)(read_otides.o) \
	$(ARC_LIB)(revrec.o) \
	$(ARC_LIB)(satprop.o) \
	$(ARC_LIB)(sbfn.o) \
	$(ARC_LIB)(sbfn1.o) \
	$(ARC_LIB)(sbout.o) \
	$(ARC_LIB)(sclcof.o) \
	$(ARC_LIB)(sclcof1.o) \
	$(ARC_LIB)(set_model_controls.o) \
	$(ARC_LIB)(shadow.o) \
	$(ARC_LIB)(shwprt.o) \
	$(ARC_LIB)(srpfgrid.o) \
	$(ARC_LIB)(srpfmod.o) \
	$(ARC_LIB)(start_int.o) \
	$(ARC_LIB)(SURFBOXW.o) \
	$(ARC_LIB)(UCLSRPeditIIA.o) \
	$(ARC_LIB)(wrthed.o) \
	$(ARC_LIB)(write_arcouthd.o) \
	$(ARC_LIB)(generalrel.o) \
	$(ARC_LIB)(print_csnm.o) \
	$(ARC_LIB)(antpwr_base.o)
	ranlib arc_lib.a

$(ARC_LIB)(adam.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(arcmrg.o): ../includes/arc.h ../includes/dimpar.h ../includes/units.h
$(ARC_LIB)(binoml.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(calcof.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(check_gmodels.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
$(ARC_LIB)(dtwopi.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(earthradTUM.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
$(ARC_LIB)(eclout.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(ertorb.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(eval.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(filopn.o): ../includes/arc.h ../includes/dimpar.h ../includes/units.h
$(ARC_LIB)(get_sat_info.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(hofrad.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(igs92.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(egm96.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(egm08.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(iers_etides.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(init.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(keplr.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(read_arc_batch.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(read_otides.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(revrec.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(sbfn.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
$(ARC_LIB)(sbfn1.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(sbout.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(sclcof.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(set_model_controls.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(shadow.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(shwprt.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(srpfgrid.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(start_int.o): ../includes/arc.h ../includes/dimpar.h
$(ARC_LIB)(wrthed.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h
$(ARC_LIB)(write_arcouthd.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h

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
