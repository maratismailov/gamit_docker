# Makefile for glinit on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:32:30 UTC 2021

all: glinit



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

GLINIT_LIB = glinit_lib.a
SLIB = ../globk/globk_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glinit : glinit_main.f $(GLINIT_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) glinit_main.f $(GLINIT_LIB) $(SLIB) -o glinit
	rm -f glinit_main.o

$(GLINIT_LIB) : \
	$(GLINIT_LIB)(glinit.o) \
	$(GLINIT_LIB)(decode_glini.o) \
	$(GLINIT_LIB)(save_epoch.o) \
	$(GLINIT_LIB)(save_apr_val.o) \
	$(GLINIT_LIB)(get_aprioris.o) \
	$(GLINIT_LIB)(save_apr_inf.o) \
	$(GLINIT_LIB)(check_apr.o) \
	$(GLINIT_LIB)(sort_epochs.o) \
	$(GLINIT_LIB)(sort_apriori.o) \
	$(GLINIT_LIB)(sort_eq.o) \
	$(GLINIT_LIB)(write_sort.o) \
	$(GLINIT_LIB)(global_init.o) \
	$(GLINIT_LIB)(solution_inf.o) \
	$(GLINIT_LIB)(emasize.o) \
	$(GLINIT_LIB)(read_eq_file.o) \
	$(GLINIT_LIB)(upd_svs_file.o) \
	$(GLINIT_LIB)(glb_glnt_mem.o)
	ranlib glinit_lib.a

$(GLINIT_LIB)(glinit.o): ../includes/glb_hdr_def.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(decode_glini.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(save_epoch.o): ../includes/kalman_param.h
$(GLINIT_LIB)(save_apr_val.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(get_aprioris.o): ../glinit/qsvi_rec.h ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/svinf_def.h
$(GLINIT_LIB)(save_apr_inf.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(sort_epochs.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(sort_apriori.o): ../includes/const_param.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(sort_eq.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(write_sort.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(global_init.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(solution_inf.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h
$(GLINIT_LIB)(read_eq_file.o): ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(upd_svs_file.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLINIT_LIB)(glb_glnt_mem.o): ../includes/kalman_param.h

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
