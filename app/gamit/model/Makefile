# Makefile for model on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:21:33 UTC 2021

all: model



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

MODEL_LIB = model_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a ../lib/libsofa.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

model : model.f $(MODEL_LIB) $(SLIB) ../../libraries/includes/const_param.h ../includes/arc.h ../includes/dimpar.h ../includes/errflg.h ../includes/global.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) model.f $(MODEL_LIB) $(SLIB) -o model
	rm -f model.o

$(MODEL_LIB) : \
	$(MODEL_LIB)(get_atmdel.o) \
	$(MODEL_LIB)(avclck.o) \
	$(MODEL_LIB)(az_elev.o) \
	$(MODEL_LIB)(call_mag.o) \
	$(MODEL_LIB)(callppt.o) \
	$(MODEL_LIB)(cfout.o) \
	$(MODEL_LIB)(chdred.o) \
	$(MODEL_LIB)(check_met.o) \
	$(MODEL_LIB)(corsim.o) \
	$(MODEL_LIB)(cortyp.o) \
	$(MODEL_LIB)(dchao.o) \
	$(MODEL_LIB)(dipole_comp.o) \
	$(MODEL_LIB)(eopart.o) \
	$(MODEL_LIB)(epoch_clk.o) \
	$(MODEL_LIB)(etide.o) \
	$(MODEL_LIB)(gdetic.o) \
	$(MODEL_LIB)(get_antinfo.o) \
	$(MODEL_LIB)(get_antpcv.o) \
	$(MODEL_LIB)(get_svantpcv.o) \
	$(MODEL_LIB)(get_dcb.o) \
	$(MODEL_LIB)(get_dcb2.o) \
	$(MODEL_LIB)(get_met_source.o) \
	$(MODEL_LIB)(svbody_coords.o) \
	$(MODEL_LIB)(gmf.o) \
	$(MODEL_LIB)(iono_corr.o) \
	$(MODEL_LIB)(iondel.o) \
	$(MODEL_LIB)(height.o) \
	$(MODEL_LIB)(iers2003_etide.o) \
	$(MODEL_LIB)(iers2010_etide.o) \
	$(MODEL_LIB)(imfh1p0.o) \
	$(MODEL_LIB)(imfw2.o) \
	$(MODEL_LIB)(lininterp.o) \
	$(MODEL_LIB)(nmfw2.o) \
	$(MODEL_LIB)(obsmod.o) \
	$(MODEL_LIB)(obsred.o) \
	$(MODEL_LIB)(open.o) \
	$(MODEL_LIB)(mversn.o) \
	$(MODEL_LIB)(nmfh2p1.o) \
	$(MODEL_LIB)(mag.o) \
	$(MODEL_LIB)(mag_dipole.o) \
	$(MODEL_LIB)(mag11.o) \
	$(MODEL_LIB)(mag12.o) \
	$(MODEL_LIB)(mag13.o) \
	$(MODEL_LIB)(igrf_util.o) \
	$(MODEL_LIB)(partl.o) \
	$(MODEL_LIB)(phasecc.o) \
	$(MODEL_LIB)(pierce_pt.o) \
	$(MODEL_LIB)(poly_clk.o) \
	$(MODEL_LIB)(ptide.o) \
	$(MODEL_LIB)(read_batch.o) \
	$(MODEL_LIB)(readi.o) \
	$(MODEL_LIB)(read_ionex.o) \
	$(MODEL_LIB)(readu.o) \
	$(MODEL_LIB)(read_metrnx.o) \
	$(MODEL_LIB)(setup.o) \
	$(MODEL_LIB)(shftor.o) \
	$(MODEL_LIB)(simred.o) \
	$(MODEL_LIB)(sitcor.o) \
	$(MODEL_LIB)(slant_tec.o) \
	$(MODEL_LIB)(svant.o) \
	$(MODEL_LIB)(update_coords.o) \
	$(MODEL_LIB)(vec_xyz.o) \
	$(MODEL_LIB)(vmf1.o) \
	$(MODEL_LIB)(wchao.o) \
	$(MODEL_LIB)(whdred.o) \
	$(MODEL_LIB)(wpress.o) \
	$(MODEL_LIB)(wrthed.o)
	ranlib model_lib.a

$(MODEL_LIB)(get_atmdel.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(avclck.o): ../includes/dimpar.h ../includes/errflg.h
$(MODEL_LIB)(az_elev.o): ../includes/dimpar.h
$(MODEL_LIB)(cfout.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(chdred.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(dipole_comp.o): ../includes/dimpar.h
$(MODEL_LIB)(eopart.o): ../includes/dimpar.h
$(MODEL_LIB)(epoch_clk.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/errflg.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(etide.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(get_antinfo.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(get_antpcv.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(get_svantpcv.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h
$(MODEL_LIB)(get_dcb.o): ../includes/dimpar.h
$(MODEL_LIB)(get_met_source.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(svbody_coords.o): ../includes/dimpar.h
$(MODEL_LIB)(iondel.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(obsmod.o): ../includes/dimpar.h
$(MODEL_LIB)(obsred.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(open.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(partl.o): ../includes/dimpar.h
$(MODEL_LIB)(phasecc.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(poly_clk.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(read_batch.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(readi.o): ../includes/dimpar.h
$(MODEL_LIB)(read_ionex.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(readu.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(read_metrnx.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(setup.o): ../../libraries/includes/freq_def.h ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(simred.o): ../includes/dimpar.h ../includes/model.h
$(MODEL_LIB)(sitcor.o): ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(svant.o): ../includes/dimpar.h
$(MODEL_LIB)(update_coords.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(MODEL_LIB)(wrthed.o): ../../libraries/includes/const_param.h ../includes/dimpar.h ../includes/global.h ../includes/model.h ../includes/units.h

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
