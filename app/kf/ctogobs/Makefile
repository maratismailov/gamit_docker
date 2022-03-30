# Makefile for ctogobs on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:35:24 UTC 2021

all: ctogobs_lib autcln ctogobs dd_sd sd_ow diff_igs merge_igs_clk



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

CTOGOBS_LIB = ctogobs_lib.a
SLIB = ./ctogobs_lib.a ../Ghandlers/cf_hndl_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

ctogobs_lib : $(CTOGOBS_LIB)

autcln : ctogobs.f  $(SLIB) ../includes/cfile_def.h ../includes/gobs_def.h ../includes/kalman_param.h ctog_cmds_bd.f ctogobs_com.h
	$(FTN) $(FFLAGS) ctogobs.f  $(SLIB) -o autcln
	rm -f ctogobs.o

ctogobs : ctogobs.f  $(SLIB) ../includes/cfile_def.h ../includes/gobs_def.h ../includes/kalman_param.h ctog_cmds_bd.f ctogobs_com.h
	$(FTN) $(FFLAGS) ctogobs.f  $(SLIB) -o ctogobs
	rm -f ctogobs.o

dd_sd : dd_sd.f  $(SLIB) 
	$(FTN) $(FFLAGS) dd_sd.f  $(SLIB) -o dd_sd
	rm -f dd_sd.o

sd_ow : sd_ow.f  $(SLIB) 
	$(FTN) $(FFLAGS) sd_ow.f  $(SLIB) -o sd_ow
	rm -f sd_ow.o

diff_igs : diff_igs.f  $(SLIB) 
	$(FTN) $(FFLAGS) diff_igs.f  $(SLIB) -o diff_igs
	rm -f diff_igs.o

merge_igs_clk : merge_igs_clk.f  $(SLIB) ../includes/kalman_param.h merge_igs_clk.h
	$(FTN) $(FFLAGS) merge_igs_clk.f  $(SLIB) -o merge_igs_clk
	rm -f merge_igs_clk.o

$(CTOGOBS_LIB) : \
	$(CTOGOBS_LIB)(ctog_utils.o) \
	$(CTOGOBS_LIB)(ctog_utm.o) \
	$(CTOGOBS_LIB)(ctog_ut2.o) \
	$(CTOGOBS_LIB)(upd_cf_sd.o) \
	$(CTOGOBS_LIB)(flag_np.o) \
	$(CTOGOBS_LIB)(trim_oneways.o) \
	$(CTOGOBS_LIB)(ctogobs_dd.o) \
	$(CTOGOBS_LIB)(ctogobs_dd2.o) \
	$(CTOGOBS_LIB)(make_sd_file.o) \
	$(CTOGOBS_LIB)(prescan_clk.o) \
	$(CTOGOBS_LIB)(set_pre_edit.o) \
	$(CTOGOBS_LIB)(chk_scan_ed.o) \
	$(CTOGOBS_LIB)(init_scan_ed.o) \
	$(CTOGOBS_LIB)(rd_ctog_cmds.o) \
	$(CTOGOBS_LIB)(prc_ctog_cmd.o) \
	$(CTOGOBS_LIB)(set_dd_ret.o) \
	$(CTOGOBS_LIB)(wr_dd_scan.o) \
	$(CTOGOBS_LIB)(rep_params.o) \
	$(CTOGOBS_LIB)(ctogobs_scan.o) \
	$(CTOGOBS_LIB)(set_scan_lis.o) \
	$(CTOGOBS_LIB)(flag_gaps.o) \
	$(CTOGOBS_LIB)(remove_fb.o) \
	$(CTOGOBS_LIB)(push_ddbf.o) \
	$(CTOGOBS_LIB)(report_bf.o) \
	$(CTOGOBS_LIB)(find_dd_slip.o) \
	$(CTOGOBS_LIB)(get_caprog.o) \
	$(CTOGOBS_LIB)(read_dfile.o) \
	$(CTOGOBS_LIB)(save_toofew.o) \
	$(CTOGOBS_LIB)(force_lc_pol.o) \
	$(CTOGOBS_LIB)(gen_gobs.o) \
	$(CTOGOBS_LIB)(comp_dphs.o) \
	$(CTOGOBS_LIB)(proc_phsfin.o) \
	$(CTOGOBS_LIB)(ctog_ow.o) \
	$(CTOGOBS_LIB)(pf_ow_sf.o) \
	$(CTOGOBS_LIB)(est_bps.o) \
	$(CTOGOBS_LIB)(flat_dd.o) \
	$(CTOGOBS_LIB)(wr_igs_clk.o) \
	$(CTOGOBS_LIB)(align_phs.o) \
	$(CTOGOBS_LIB)(fit_igs_clk.o) \
	$(CTOGOBS_LIB)(zero_dd_wl.o) \
	$(CTOGOBS_LIB)(est_wlfull.o) \
	$(CTOGOBS_LIB)(est_dd_wl.o) \
	$(CTOGOBS_LIB)(set_freqs.o) \
	$(CTOGOBS_LIB)(comp_ion.o) \
	$(CTOGOBS_LIB)(fit_glonass_clks.o) \
	$(CTOGOBS_LIB)(prescan_msj.o) \
	$(CTOGOBS_LIB)(verify_scan_dd.o) \
	$(CTOGOBS_LIB)(trim_shortseg.o)
	ranlib ctogobs_lib.a

$(CTOGOBS_LIB)(ctog_utils.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/gobs_def.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(ctog_utm.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/gobs_header.h ../includes/kalman_param.h ../includes/mfile_def.h ctogobs_com.h
$(CTOGOBS_LIB)(ctog_ut2.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(upd_cf_sd.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(flag_np.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(trim_oneways.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(ctogobs_dd.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(ctogobs_dd2.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(make_sd_file.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(prescan_clk.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(set_pre_edit.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(chk_scan_ed.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(rd_ctog_cmds.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(prc_ctog_cmd.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(rep_params.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(ctogobs_scan.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(set_scan_lis.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(flag_gaps.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(remove_fb.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(push_ddbf.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(report_bf.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(find_dd_slip.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(force_lc_pol.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(gen_gobs.o): ../includes/cfile_def.h ../includes/gobs_def.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(comp_dphs.o): ../includes/cfile_def.h ../includes/kalman_param.h ../includes/mfile_def.h ctogobs_com.h
$(CTOGOBS_LIB)(proc_phsfin.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(ctog_ow.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(pf_ow_sf.o): ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(est_bps.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(flat_dd.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(wr_igs_clk.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ../includes/mfile_def.h ctogobs_com.h
$(CTOGOBS_LIB)(align_phs.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(fit_igs_clk.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ../includes/mfile_def.h ctogobs_com.h
$(CTOGOBS_LIB)(zero_dd_wl.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(est_wlfull.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(est_dd_wl.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(set_freqs.o): ../includes/cfile_def.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(comp_ion.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(fit_glonass_clks.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(prescan_msj.o): ../includes/cfile_def.h ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(verify_scan_dd.o): ../includes/const_param.h ../includes/kalman_param.h ctogobs_com.h
$(CTOGOBS_LIB)(trim_shortseg.o): ../includes/kalman_param.h ctogobs_com.h

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
