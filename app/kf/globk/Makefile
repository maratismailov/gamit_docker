# Makefile for globk on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:30:49 UTC 2021

all: globs



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

GLOBK_LIB = globk_lib.a
SLIB = ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

globs : globk.f $(GLOBK_LIB) $(SLIB) ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h globk_cmd_bd.f
	$(FTN) $(FFLAGS) globk.f $(GLOBK_LIB) $(SLIB) -o globs
	rm -f globk.o

$(GLOBK_LIB) : \
	$(GLOBK_LIB)(decode_globk.o) \
	$(GLOBK_LIB)(decode_prt.o) \
	$(GLOBK_LIB)(proper_glb_r.o) \
	$(GLOBK_LIB)(read_glb_mar.o) \
	$(GLOBK_LIB)(process_glb_.o) \
	$(GLOBK_LIB)(assign_ema_v.o) \
	$(GLOBK_LIB)(run_glinit.o) \
	$(GLOBK_LIB)(run_glfor.o) \
	$(GLOBK_LIB)(run_glsave.o) \
	$(GLOBK_LIB)(run_glbak.o) \
	$(GLOBK_LIB)(run_glout.o) \
	$(GLOBK_LIB)(run_glorg.o) \
	$(GLOBK_LIB)(gen_org_name.o) \
	$(GLOBK_LIB)(get_apr_posi.o) \
	$(GLOBK_LIB)(try_upd_apr.o) \
	$(GLOBK_LIB)(glb_param_li.o) \
	$(GLOBK_LIB)(clear_glb_pa.o) \
	$(GLOBK_LIB)(glb_par_coun.o) \
	$(GLOBK_LIB)(check_trans.o) \
	$(GLOBK_LIB)(check_glb_ma.o) \
	$(GLOBK_LIB)(globk_header.o) \
	$(GLOBK_LIB)(glfor_mem.o) \
	$(GLOBK_LIB)(glb_for_map.o) \
	$(GLOBK_LIB)(glb_bak_map.o) \
	$(GLOBK_LIB)(globk_extras.o) \
	$(GLOBK_LIB)(decode_nonse.o) \
	$(GLOBK_LIB)(add_GGV.o)
	ranlib globk_lib.a

$(GLOBK_LIB)(decode_globk.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(read_glb_mar.o): ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(process_glb_.o): ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glinit.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glfor.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glsave.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glbak.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glout.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(run_glorg.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(get_apr_posi.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(try_upd_apr.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(glb_param_li.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(clear_glb_pa.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(glb_par_coun.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(check_trans.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(glfor_mem.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(glb_for_map.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLOBK_LIB)(glb_bak_map.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLOBK_LIB)(globk_extras.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLOBK_LIB)(decode_nonse.o): ../includes/globk_common.h ../includes/kalman_param.h

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
