# Makefile for glbak on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:33:48 UTC 2021

all: glbak



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

GLBAK_LIB = glbak_lib.a
SLIB = ../glorg/glorg_lib.a ../glout/glout_lib.a ../glfor/glfor_lib.a ../glsave/glsave_lib.a ../globk/globk_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glbak : glbak_main.f $(GLBAK_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) glbak_main.f $(GLBAK_LIB) $(SLIB) -o glbak
	rm -f glbak_main.o

$(GLBAK_LIB) : \
	$(GLBAK_LIB)(glbak.o) \
	$(GLBAK_LIB)(glb_bak_filt.o) \
	$(GLBAK_LIB)(glb_bak_writ.o) \
	$(GLBAK_LIB)(set_mar_prt.o) \
	$(GLBAK_LIB)(glb_average.o) \
	$(GLBAK_LIB)(init_bak_cov.o) \
	$(GLBAK_LIB)(glb_residual.o) \
	$(GLBAK_LIB)(write_bak_he.o) \
	$(GLBAK_LIB)(write_loc_pa.o) \
	$(GLBAK_LIB)(out_glbl.o) \
	$(GLBAK_LIB)(check_est.o) \
	$(GLBAK_LIB)(clear_loc_pa.o) \
	$(GLBAK_LIB)(get_bak_apr.o) \
	$(GLBAK_LIB)(parm_to_obs.o) \
	$(GLBAK_LIB)(correct_obs.o) \
	$(GLBAK_LIB)(glb_inc_post.o) \
	$(GLBAK_LIB)(inc_bak_stat.o) \
	$(GLBAK_LIB)(keep_loc.o) \
	$(GLBAK_LIB)(glb_inc_bak.o) \
	$(GLBAK_LIB)(read_glbak_orgcmd.o) \
	$(GLBAK_LIB)(decode_glbak_run.o)
	ranlib glbak_lib.a

$(GLBAK_LIB)(glbak.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_bak_filt.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_bak_writ.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(set_mar_prt.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_average.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(init_bak_cov.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_residual.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLBAK_LIB)(write_bak_he.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(write_loc_pa.o): ../includes/const_param.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(get_bak_apr.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_inc_post.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLBAK_LIB)(inc_bak_stat.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLBAK_LIB)(glb_inc_bak.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h
$(GLBAK_LIB)(read_glbak_orgcmd.o): ../includes/globk_common.h ../includes/glorg_common.h ../includes/kalman_param.h

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
