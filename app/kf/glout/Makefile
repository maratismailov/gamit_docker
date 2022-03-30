# Makefile for glout on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:32:00 UTC 2021

all: glout



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

GLOUT_LIB = glout_lib.a
SLIB = ../glfor/glfor_lib.a ../globk/globk_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glout : glout_main.f $(GLOUT_LIB) $(SLIB) ../globk/globk_cmd_bd.f
	$(FTN) $(FFLAGS) glout_main.f $(GLOUT_LIB) $(SLIB) -o glout
	rm -f glout_main.o

$(GLOUT_LIB) : \
	$(GLOUT_LIB)(glout.o) \
	$(GLOUT_LIB)(glouo.o) \
	$(GLOUT_LIB)(glous.o) \
	$(GLOUT_LIB)(glouw.o) \
	$(GLOUT_LIB)(decode_glout.o) \
	$(GLOUT_LIB)(glb_out_map.o) \
	$(GLOUT_LIB)(write_glb_he.o) \
	$(GLOUT_LIB)(write_glb_pa.o) \
	$(GLOUT_LIB)(write_glb_ba.o) \
	$(GLOUT_LIB)(write_glb_bc.o) \
	$(GLOUT_LIB)(write_glb_su.o) \
	$(GLOUT_LIB)(write_glb_sp.o) \
	$(GLOUT_LIB)(summ_eq.o) \
	$(GLOUT_LIB)(write_iers.o) \
	$(GLOUT_LIB)(out_glbp.o) \
	$(GLOUT_LIB)(out_type_sta.o) \
	$(GLOUT_LIB)(get_from_ema.o) \
	$(GLOUT_LIB)(put_to_ema.o) \
	$(GLOUT_LIB)(write_glb_co.o) \
	$(GLOUT_LIB)(comp_nut_amp.o) \
	$(GLOUT_LIB)(get_nut_coef.o) \
	$(GLOUT_LIB)(get_sd_coeff.o) \
	$(GLOUT_LIB)(report_eq.o) \
	$(GLOUT_LIB)(report_noapr.o) \
	$(GLOUT_LIB)(qoutbas.o) \
	$(GLOUT_LIB)(elem.o) \
	$(GLOUT_LIB)(report_nonse.o)
	ranlib glout_lib.a

$(GLOUT_LIB)(glout.o): ../includes/globk_common.h ../includes/kalman_param.h ../includes/sd_common.h
$(GLOUT_LIB)(glouo.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(glous.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(glouw.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(glb_out_map.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_he.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_pa.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/sd_common.h
$(GLOUT_LIB)(write_glb_ba.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_bc.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_su.o): ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_sp.o): ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(summ_eq.o): ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h ../includes/sd_common.h
$(GLOUT_LIB)(write_iers.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(out_glbp.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLOUT_LIB)(out_type_sta.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(write_glb_co.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(get_nut_coef.o): ../includes/nut_eval.h
$(GLOUT_LIB)(report_eq.o): ../includes/const_param.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(report_noapr.o): ../includes/const_param.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(qoutbas.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOUT_LIB)(report_nonse.o): ../includes/globk_common.h ../includes/kalman_param.h

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
