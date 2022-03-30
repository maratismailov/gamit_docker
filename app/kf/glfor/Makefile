# Makefile for glfor on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:31:17 UTC 2021

all: glfor



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

GLFOR_LIB = glfor_lib.a
SLIB = ../globk/globk_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

glfor : glfor_main.f $(GLFOR_LIB) $(SLIB) ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
	$(FTN) $(FFLAGS) glfor_main.f $(GLFOR_LIB) $(SLIB) -o glfor
	rm -f glfor_main.o

$(GLFOR_LIB) : \
	$(GLFOR_LIB)(glfor.o) \
	$(GLFOR_LIB)(decode_glfor.o) \
	$(GLFOR_LIB)(get_param_co.o) \
	$(GLFOR_LIB)(glb_o_minus_.o) \
	$(GLFOR_LIB)(compute_OC.o) \
	$(GLFOR_LIB)(site_OC.o) \
	$(GLFOR_LIB)(source_OC.o) \
	$(GLFOR_LIB)(glb_partials.o) \
	$(GLFOR_LIB)(clear_pnt.o) \
	$(GLFOR_LIB)(compute_part.o) \
	$(GLFOR_LIB)(site_partial.o) \
	$(GLFOR_LIB)(site_rate_pa.o) \
	$(GLFOR_LIB)(pmu_part.o) \
	$(GLFOR_LIB)(add_partial.o) \
	$(GLFOR_LIB)(source_parti.o) \
	$(GLFOR_LIB)(nut_ang_part.o) \
	$(GLFOR_LIB)(scan_parm.o) \
	$(GLFOR_LIB)(nut_coeff_pa.o) \
	$(GLFOR_LIB)(gamp_parts.o) \
	$(GLFOR_LIB)(get_max_glbd.o) \
	$(GLFOR_LIB)(coeff_OC.o) \
	$(GLFOR_LIB)(glb_for_filt.o) \
	$(GLFOR_LIB)(glb_predict.o) \
	$(GLFOR_LIB)(out_cov.o) \
	$(GLFOR_LIB)(get_trans_co.o) \
	$(GLFOR_LIB)(init_glb_cov.o) \
	$(GLFOR_LIB)(glb_kalman_g.o) \
	$(GLFOR_LIB)(glb_dot_part.o) \
	$(GLFOR_LIB)(glb_update.o) \
	$(GLFOR_LIB)(eq_glb_pred.o) \
	$(GLFOR_LIB)(eq_effected.o) \
	$(GLFOR_LIB)(find_preq.o) \
	$(GLFOR_LIB)(glb_inc_pref.o) \
	$(GLFOR_LIB)(remove_param.o) \
	$(GLFOR_LIB)(codes_to_par.o) \
	$(GLFOR_LIB)(glb_upd_apr.o) \
	$(GLFOR_LIB)(rn_app_dpos.o) \
	$(GLFOR_LIB)(get_mar_svs.o) \
	$(GLFOR_LIB)(read_pmu_inp.o) \
	$(GLFOR_LIB)(print_pmu.o) \
	$(GLFOR_LIB)(app_ptide.o) \
	$(GLFOR_LIB)(mul_pmu_pred.o) \
	$(GLFOR_LIB)(eval_nonsec.o) \
	$(GLFOR_LIB)(change_load.o) \
	$(GLFOR_LIB)(ptd_status.o)
	ranlib glfor_lib.a

$(GLFOR_LIB)(glfor.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(decode_glfor.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLFOR_LIB)(get_param_co.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_o_minus_.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(compute_OC.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h ../includes/sd_common.h
$(GLFOR_LIB)(site_OC.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(source_OC.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_partials.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(clear_pnt.o): ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLFOR_LIB)(compute_part.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(site_partial.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(site_rate_pa.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(pmu_part.o): ../includes/const_param.h
$(GLFOR_LIB)(source_parti.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(nut_ang_part.o): ../includes/const_param.h
$(GLFOR_LIB)(scan_parm.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLFOR_LIB)(get_max_glbd.o): ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_for_filt.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_predict.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(init_glb_cov.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_kalman_g.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_update.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/kalman_param.h
$(GLFOR_LIB)(eq_glb_pred.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(eq_effected.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(find_preq.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_inc_pref.o): ../includes/glb_hdr_def.h ../includes/globk_cntl.h ../includes/globk_markov.h ../includes/kalman_param.h ../includes/sln_def.h
$(GLFOR_LIB)(remove_param.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(glb_upd_apr.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(rn_app_dpos.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(get_mar_svs.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(read_pmu_inp.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(print_pmu.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(app_ptide.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(mul_pmu_pred.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h
$(GLFOR_LIB)(eval_nonsec.o): ../includes/const_param.h
$(GLFOR_LIB)(change_load.o): ../includes/glb_hdr_def.h ../includes/globk_common.h ../includes/kalman_param.h

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
