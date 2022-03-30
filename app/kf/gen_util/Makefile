# Makefile for gen_util on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:28:35 UTC 2021

all: gen_util



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

GEN_UTIL_LIB = gen_util_lib.a
SLIB = ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

gen_util : $(GEN_UTIL_LIB)

$(GEN_UTIL_LIB) : \
	$(GEN_UTIL_LIB)(gstut.o) \
	$(GEN_UTIL_LIB)(gst_jd.o) \
	$(GEN_UTIL_LIB)(epoc_8.o) \
	$(GEN_UTIL_LIB)(smnv8.o) \
	$(GEN_UTIL_LIB)(wt_filter_fw.o) \
	$(GEN_UTIL_LIB)(wt_gauss_fil.o) \
	$(GEN_UTIL_LIB)(lagrange_int.o) \
	$(GEN_UTIL_LIB)(short_ut1.o) \
	$(GEN_UTIL_LIB)(short_lod.o) \
	$(GEN_UTIL_LIB)(fjldy_8.o) \
	$(GEN_UTIL_LIB)(baseline_num.o) \
	$(GEN_UTIL_LIB)(gauss_wrms.o) \
	$(GEN_UTIL_LIB)(int_to_char.o) \
	$(GEN_UTIL_LIB)(jel.o) \
	$(GEN_UTIL_LIB)(mas_to_dms.o) \
	$(GEN_UTIL_LIB)(sat_press.o) \
	$(GEN_UTIL_LIB)(compress_rc.o) \
	$(GEN_UTIL_LIB)(cross_prod.o) \
	$(GEN_UTIL_LIB)(clean_file.o) \
	$(GEN_UTIL_LIB)(rotate_geod.o) \
	$(GEN_UTIL_LIB)(loc_to_geod.o) \
	$(GEN_UTIL_LIB)(dot3.o) \
	$(GEN_UTIL_LIB)(diagonalize.o) \
	$(GEN_UTIL_LIB)(check_covar.o) \
	$(GEN_UTIL_LIB)(rotat_crd.o) \
	$(GEN_UTIL_LIB)(report_break.o) \
	$(GEN_UTIL_LIB)(switch.o) \
	$(GEN_UTIL_LIB)(strip_name.o) \
	$(GEN_UTIL_LIB)(out_run_time.o) \
	$(GEN_UTIL_LIB)(base_entnum.o) \
	$(GEN_UTIL_LIB)(bak_file_map.o) \
	$(GEN_UTIL_LIB)(base_data.o) \
	$(GEN_UTIL_LIB)(num_step.o) \
	$(GEN_UTIL_LIB)(cal_date.o) \
	$(GEN_UTIL_LIB)(clear_close.o) \
	$(GEN_UTIL_LIB)(close_cont.o) \
	$(GEN_UTIL_LIB)(clear_norm.o) \
	$(GEN_UTIL_LIB)(increment_de.o) \
	$(GEN_UTIL_LIB)(decode_mar.o) \
	$(GEN_UTIL_LIB)(decode_optio.o) \
	$(GEN_UTIL_LIB)(set_use.o) \
	$(GEN_UTIL_LIB)(con_mas.o) \
	$(GEN_UTIL_LIB)(duse.o) \
	$(GEN_UTIL_LIB)(construct_th.o) \
	$(GEN_UTIL_LIB)(add_cont.o) \
	$(GEN_UTIL_LIB)(add_struct.o) \
	$(GEN_UTIL_LIB)(add_site_con.o) \
	$(GEN_UTIL_LIB)(add_ion.o) \
	$(GEN_UTIL_LIB)(add_medium.o) \
	$(GEN_UTIL_LIB)(apply_rot.o) \
	$(GEN_UTIL_LIB)(apply_trans.o) \
	$(GEN_UTIL_LIB)(add_clk_brk.o) \
	$(GEN_UTIL_LIB)(etd_partial.o) \
	$(GEN_UTIL_LIB)(eor_partial.o) \
	$(GEN_UTIL_LIB)(plm.o) \
	$(GEN_UTIL_LIB)(plmdz.o) \
	$(GEN_UTIL_LIB)(atm_az_parti.o) \
	$(GEN_UTIL_LIB)(atm_map_part.o) \
	$(GEN_UTIL_LIB)(dry_saas.o) \
	$(GEN_UTIL_LIB)(dry_hgt.o) \
	$(GEN_UTIL_LIB)(dry_hgt_site.o) \
	$(GEN_UTIL_LIB)(dry_const.o) \
	$(GEN_UTIL_LIB)(wet_saas.o) \
	$(GEN_UTIL_LIB)(Wet_WVR_db.o) \
	$(GEN_UTIL_LIB)(wet_const.o) \
	$(GEN_UTIL_LIB)(dry_mari.o) \
	$(GEN_UTIL_LIB)(dry_cfa.o) \
	$(GEN_UTIL_LIB)(dry_mit.o) \
	$(GEN_UTIL_LIB)(dry_chao.o) \
	$(GEN_UTIL_LIB)(dry_mtt.o) \
	$(GEN_UTIL_LIB)(wet_chao.o) \
	$(GEN_UTIL_LIB)(wet_csc.o) \
	$(GEN_UTIL_LIB)(wet_mari.o) \
	$(GEN_UTIL_LIB)(wet_cfa.o) \
	$(GEN_UTIL_LIB)(wet_mit.o) \
	$(GEN_UTIL_LIB)(wet_mtt.o) \
	$(GEN_UTIL_LIB)(met_seasonal.o) \
	$(GEN_UTIL_LIB)(initialize_a.o) \
	$(GEN_UTIL_LIB)(wild_card.o) \
	$(GEN_UTIL_LIB)(wild_date.o) \
	$(GEN_UTIL_LIB)(lastex.o) \
	$(GEN_UTIL_LIB)(apply_edits.o) \
	$(GEN_UTIL_LIB)(reset_data_f.o) \
	$(GEN_UTIL_LIB)(wet_ifad.o) \
	$(GEN_UTIL_LIB)(wet_pressure.o) \
	$(GEN_UTIL_LIB)(phi_function.o) \
	$(GEN_UTIL_LIB)(get_gen_rel.o) \
	$(GEN_UTIL_LIB)(get_ema_4.o) \
	$(GEN_UTIL_LIB)(get_ema_8.o) \
	$(GEN_UTIL_LIB)(get_bit_stat.o) \
	$(GEN_UTIL_LIB)(var_comp.o) \
	$(GEN_UTIL_LIB)(write_line.o) \
	$(GEN_UTIL_LIB)(mov_cov.o) \
	$(GEN_UTIL_LIB)(out_loc_corr.o) \
	$(GEN_UTIL_LIB)(fix_origin.o) \
	$(GEN_UTIL_LIB)(glb_stabiliz.o) \
	$(GEN_UTIL_LIB)(ltog_map.o) \
	$(GEN_UTIL_LIB)(gtol_map.o) \
	$(GEN_UTIL_LIB)(fix_system.o) \
	$(GEN_UTIL_LIB)(comp_hchi.o) \
	$(GEN_UTIL_LIB)(tran_to_org.o) \
	$(GEN_UTIL_LIB)(pmu_main_par.o) \
	$(GEN_UTIL_LIB)(sd_util.o) \
	$(GEN_UTIL_LIB)(map_part.o) \
	$(GEN_UTIL_LIB)(fund_arg.o) \
	$(GEN_UTIL_LIB)(planetary_ar.o) \
	$(GEN_UTIL_LIB)(read_spc.o) \
	$(GEN_UTIL_LIB)(compute_std.o) \
	$(GEN_UTIL_LIB)(read_nut_ser.o) \
	$(GEN_UTIL_LIB)(list_markov.o) \
	$(GEN_UTIL_LIB)(ocean_util.o) \
	$(GEN_UTIL_LIB)(set_oce_avai.o) \
	$(GEN_UTIL_LIB)(force_parms.o) \
	$(GEN_UTIL_LIB)(esort.o) \
	$(GEN_UTIL_LIB)(write_pmu_co.o) \
	$(GEN_UTIL_LIB)(apply_cond.o) \
	$(GEN_UTIL_LIB)(eval_dist.o) \
	$(GEN_UTIL_LIB)(frame_to_fra.o) \
	$(GEN_UTIL_LIB)(svel_to_code.o) \
	$(GEN_UTIL_LIB)(get_institut.o) \
	$(GEN_UTIL_LIB)(hsver.o) \
	$(GEN_UTIL_LIB)(precess.o) \
	$(GEN_UTIL_LIB)(mdian2.o) \
	$(GEN_UTIL_LIB)(comp_ptide.o) \
	$(GEN_UTIL_LIB)(gamit_mod_n.o) \
	$(GEN_UTIL_LIB)(get_leapsec.o) \
	$(GEN_UTIL_LIB)(real_stats.o) \
	$(GEN_UTIL_LIB)(datum_def.o) \
	$(GEN_UTIL_LIB)(GEOD_to_UTM.o) \
	$(GEN_UTIL_LIB)(GEOD_to_GEOD.o) \
	$(GEN_UTIL_LIB)(cov_xyz_neu.o) \
	$(GEN_UTIL_LIB)(decode_comopt.o) \
	$(GEN_UTIL_LIB)(nonsec_util.o) \
	$(GEN_UTIL_LIB)(sp3_lib.o) \
	$(GEN_UTIL_LIB)(con_util.o) \
	$(GEN_UTIL_LIB)(svcode_to_code.o)
	ranlib gen_util_lib.a

$(GEN_UTIL_LIB)(gst_jd.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(loc_to_geod.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(report_break.o): ../includes/kalman_param.h ../includes/obs_names.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(clear_norm.o): ../includes/kalman_param.h
$(GEN_UTIL_LIB)(construct_th.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(add_ion.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(add_medium.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(etd_partial.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(eor_partial.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(atm_map_part.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_saas.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_hgt.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_hgt_site.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_const.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_saas.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(Wet_WVR_db.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_const.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_mari.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_cfa.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_mit.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(dry_chao.o): ../includes/kalman_param.h ../includes/obs_data.h
$(GEN_UTIL_LIB)(dry_mtt.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_chao.o): ../includes/kalman_param.h ../includes/obs_data.h
$(GEN_UTIL_LIB)(wet_csc.o): ../includes/kalman_param.h ../includes/obs_data.h
$(GEN_UTIL_LIB)(wet_mari.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_cfa.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_mit.o): ../includes/kalman_param.h ../includes/obs_data.h
$(GEN_UTIL_LIB)(wet_mtt.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(met_seasonal.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(apply_edits.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(wet_ifad.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(get_gen_rel.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_apr.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(pmu_main_par.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(sd_util.o): ../includes/const_param.h ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_header.h ../includes/sd_common.h
$(GEN_UTIL_LIB)(map_part.o): ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_values.h
$(GEN_UTIL_LIB)(fund_arg.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(planetary_ar.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(read_spc.o): ../includes/nut_eval.h
$(GEN_UTIL_LIB)(compute_std.o): ../includes/const_param.h ../includes/nut_eval.h
$(GEN_UTIL_LIB)(read_nut_ser.o): ../includes/nut_eval.h
$(GEN_UTIL_LIB)(ocean_util.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(frame_to_fra.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(precess.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(comp_ptide.o): ../includes/const_param.h
$(GEN_UTIL_LIB)(datum_def.o): ../includes/const_param.h ../includes/utmut.h
$(GEN_UTIL_LIB)(GEOD_to_UTM.o): ../includes/const_param.h ../includes/utmut.h
$(GEN_UTIL_LIB)(GEOD_to_GEOD.o): ../includes/const_param.h ../includes/utmut.h
$(GEN_UTIL_LIB)(cov_xyz_neu.o): ../includes/globk_common.h ../includes/kalman_param.h
$(GEN_UTIL_LIB)(nonsec_util.o): ../gen_util/nonsec_util.h ../includes/const_param.h
$(GEN_UTIL_LIB)(sp3_lib.o): ../../libraries/includes/freq_def.h ../includes/const_param.h ../includes/sp3_def.h
$(GEN_UTIL_LIB)(svcode_to_code.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h

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
