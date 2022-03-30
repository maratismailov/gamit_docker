# Makefile for Khandlers on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:29:57 UTC 2021

all: KalObs_hndl_lib



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

KALOBS_HNDL_LIB = KalObs_hndl_lib.a
SLIB = 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

KalObs_hndl_lib : $(KALOBS_HNDL_LIB)

$(KALOBS_HNDL_LIB) : \
	$(KALOBS_HNDL_LIB)(out_header.o) \
	$(KALOBS_HNDL_LIB)(rw_KalObs_he.o) \
	$(KALOBS_HNDL_LIB)(rw_KalObs_bl.o) \
	$(KALOBS_HNDL_LIB)(open_KalObs.o) \
	$(KALOBS_HNDL_LIB)(close_KalObs.o) \
	$(KALOBS_HNDL_LIB)(create_KalOb.o) \
	$(KALOBS_HNDL_LIB)(corrupt_KalO.o) \
	$(KALOBS_HNDL_LIB)(report_edit.o) \
	$(KALOBS_HNDL_LIB)(wr_site_mark.o) \
	$(KALOBS_HNDL_LIB)(extended_hea.o) \
	$(KALOBS_HNDL_LIB)(open_simulat.o) \
	$(KALOBS_HNDL_LIB)(ifmar.o) \
	$(KALOBS_HNDL_LIB)(rw_solvk_blo.o) \
	$(KALOBS_HNDL_LIB)(open_solvk.o) \
	$(KALOBS_HNDL_LIB)(close_solvk.o) \
	$(KALOBS_HNDL_LIB)(create_solvk.o) \
	$(KALOBS_HNDL_LIB)(rw_common.o) \
	$(KALOBS_HNDL_LIB)(rw_covar.o) \
	$(KALOBS_HNDL_LIB)(rw_glb_covar.o) \
	$(KALOBS_HNDL_LIB)(read_kalfile.o) \
	$(KALOBS_HNDL_LIB)(rw_Globk_blo.o) \
	$(KALOBS_HNDL_LIB)(rw_Ema_block.o) \
	$(KALOBS_HNDL_LIB)(open_globk.o) \
	$(KALOBS_HNDL_LIB)(close_GLOBK.o) \
	$(KALOBS_HNDL_LIB)(create_GLOBK.o) \
	$(KALOBS_HNDL_LIB)(clear_ema.o) \
	$(KALOBS_HNDL_LIB)(rw_globk_com.o) \
	$(KALOBS_HNDL_LIB)(globk_sum.o) \
	$(KALOBS_HNDL_LIB)(rw_glb_heade.o) \
	$(KALOBS_HNDL_LIB)(rw_names_blo.o) \
	$(KALOBS_HNDL_LIB)(get_name_fro.o) \
	$(KALOBS_HNDL_LIB)(add_to_list.o) \
	$(KALOBS_HNDL_LIB)(eq_name_chng.o) \
	$(KALOBS_HNDL_LIB)(add_to_buffe.o) \
	$(KALOBS_HNDL_LIB)(decode_code.o) \
	$(KALOBS_HNDL_LIB)(switch_ema8.o) \
	$(KALOBS_HNDL_LIB)(switch_ema4.o) \
	$(KALOBS_HNDL_LIB)(save_glb_apr.o) \
	$(KALOBS_HNDL_LIB)(inc_kalver.o) \
	$(KALOBS_HNDL_LIB)(wread.o) \
	$(KALOBS_HNDL_LIB)(wwrit.o) \
	$(KALOBS_HNDL_LIB)(inc_rename.o) \
	$(KALOBS_HNDL_LIB)(read_sinf_recs.o) \
	$(KALOBS_HNDL_LIB)(read_svinf_rec.o)
	ranlib KalObs_hndl_lib.a

$(KALOBS_HNDL_LIB)(out_header.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(rw_KalObs_he.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(rw_KalObs_bl.o): ../includes/kalman_param.h ../includes/obs_values.h
$(KALOBS_HNDL_LIB)(open_KalObs.o): ../includes/kalman_param.h ../includes/obs_values.h
$(KALOBS_HNDL_LIB)(close_KalObs.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(create_KalOb.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_header.h ../includes/readin_user.h
$(KALOBS_HNDL_LIB)(report_edit.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_names.h ../includes/obs_values.h
$(KALOBS_HNDL_LIB)(wr_site_mark.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(extended_hea.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(open_simulat.o): ../includes/kalman_param.h ../includes/obs_values.h
$(KALOBS_HNDL_LIB)(rw_solvk_blo.o): ../includes/kalman_param.h ../includes/solvk_cntl.h
$(KALOBS_HNDL_LIB)(open_solvk.o): ../includes/kalman_param.h ../includes/solvk_cntl.h
$(KALOBS_HNDL_LIB)(close_solvk.o): ../includes/kalman_param.h ../includes/solvk_cntl.h
$(KALOBS_HNDL_LIB)(create_solvk.o): ../includes/kalman_param.h ../includes/solvk_common.h
$(KALOBS_HNDL_LIB)(rw_common.o): ../includes/kalman_param.h ../includes/solvk_common.h
$(KALOBS_HNDL_LIB)(rw_covar.o): ../includes/kalman_param.h ../includes/obs_values.h ../includes/solvk_cntl.h
$(KALOBS_HNDL_LIB)(rw_glb_covar.o): ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(read_kalfile.o): ../includes/forsl_comm.h ../includes/kalfile_def.h ../includes/kalman_param.h ../includes/solvk_cntl.h
$(KALOBS_HNDL_LIB)(rw_Globk_blo.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(rw_Ema_block.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(open_globk.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(close_GLOBK.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(create_GLOBK.o): ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(rw_globk_com.o): ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(globk_sum.o): ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(rw_glb_heade.o): ../includes/glb_hdr_def.h ../includes/globk_markov.h ../includes/kalman_param.h ../includes/sln_def.h
$(KALOBS_HNDL_LIB)(rw_names_blo.o): ../glinit/qsvi_rec.h ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(eq_name_chng.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(save_glb_apr.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(inc_kalver.o): ../includes/kalman_param.h ../includes/obs_header.h
$(KALOBS_HNDL_LIB)(inc_rename.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(KALOBS_HNDL_LIB)(read_sinf_recs.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/sinf_def.h
$(KALOBS_HNDL_LIB)(read_svinf_rec.o): ../glinit/qsvi_rec.h ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/svinf_def.h

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
