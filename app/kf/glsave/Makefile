# Makefile for glsave on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:33:24 UTC 2021

all: glsave



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

GLSAVE_LIB = glsave_lib.a
SLIB = ../globk/globk_lib.a ../glfor/glfor_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glsave : glsave_main.f $(GLSAVE_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) glsave_main.f $(GLSAVE_LIB) $(SLIB) -o glsave
	rm -f glsave_main.o

$(GLSAVE_LIB) : \
	$(GLSAVE_LIB)(glsave.o) \
	$(GLSAVE_LIB)(gw_glb_codes.o) \
	$(GLSAVE_LIB)(glb_apr_code.o) \
	$(GLSAVE_LIB)(glb_coeff_co.o) \
	$(GLSAVE_LIB)(glb_glb_code.o) \
	$(GLSAVE_LIB)(glb_svs_code.o) \
	$(GLSAVE_LIB)(glb_tot_solu.o) \
	$(GLSAVE_LIB)(gw_glb_heade.o) \
	$(GLSAVE_LIB)(gw_glb_names.o) \
	$(GLSAVE_LIB)(gw_glb_full.o) \
	$(GLSAVE_LIB)(gw_descripti.o) \
	$(GLSAVE_LIB)(gw_codes.o) \
	$(GLSAVE_LIB)(gw_mul_ep.o) \
	$(GLSAVE_LIB)(gw_aprioris.o) \
	$(GLSAVE_LIB)(glb_save_apr.o) \
	$(GLSAVE_LIB)(gw_soln.o) \
	$(GLSAVE_LIB)(gw_cons.o) \
	$(GLSAVE_LIB)(list_soln.o)
	ranlib glsave_lib.a

$(GLSAVE_LIB)(glsave.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h ../includes/sln_def.h
$(GLSAVE_LIB)(gw_glb_codes.o): ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h ../includes/obs_header.h
$(GLSAVE_LIB)(glb_apr_code.o): ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(glb_coeff_co.o): ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(glb_glb_code.o): ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(glb_svs_code.o): ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(glb_tot_solu.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_glb_heade.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h ../includes/obs_header.h ../includes/sln_def.h
$(GLSAVE_LIB)(gw_glb_names.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_glb_full.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_descripti.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/glsave_comm.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h ../includes/svinf_def.h
$(GLSAVE_LIB)(gw_codes.o): ../includes/glb_hdr_def.h ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_mul_ep.o): ../includes/glb_hdr_def.h ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_aprioris.o): ../includes/glb_hdr_def.h ../includes/glsave_comm.h ../includes/kalman_param.h
$(GLSAVE_LIB)(glb_save_apr.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/obs_header.h
$(GLSAVE_LIB)(gw_soln.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLSAVE_LIB)(gw_cons.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h

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
