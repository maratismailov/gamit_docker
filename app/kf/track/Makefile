# Makefile for track on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:38:43 UTC 2021

all: track_lib track



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

TRACK_LIB = track_lib.a
SLIB = ./track_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

track_lib : $(TRACK_LIB)

track : track.f  $(SLIB) ../includes/xfile_def.h track_com.h
	$(FTN) $(FFLAGS) track.f  $(SLIB) -o track
	rm -f track.o

$(TRACK_LIB) : \
	$(TRACK_LIB)(assign_con.o) \
	$(TRACK_LIB)(comp_clock.o) \
	$(TRACK_LIB)(ds_to_ymdhms.o) \
	$(TRACK_LIB)(earth_to_inert.o) \
	$(TRACK_LIB)(get_elev.o) \
	$(TRACK_LIB)(get_miss.o) \
	$(TRACK_LIB)(get_obsfile.o) \
	$(TRACK_LIB)(get_site.o) \
	$(TRACK_LIB)(geoxyz.o) \
	$(TRACK_LIB)(get_svsrun.o) \
	$(TRACK_LIB)(lowerc.o) \
	$(TRACK_LIB)(neph_to_xyz.o) \
	$(TRACK_LIB)(read_bat.o) \
	$(TRACK_LIB)(read_nav.o) \
	$(TRACK_LIB)(read_nav_clk.o) \
	$(TRACK_LIB)(ymdhms_to_ds.o) \
	$(TRACK_LIB)(save_head.o) \
	$(TRACK_LIB)(clean_sp3_clk.o) \
	$(TRACK_LIB)(init_track.o) \
	$(TRACK_LIB)(read_all_rinex.o) \
	$(TRACK_LIB)(est_pos.o) \
	$(TRACK_LIB)(theory.o) \
	$(TRACK_LIB)(search.o) \
	$(TRACK_LIB)(atm_lib.o) \
	$(TRACK_LIB)(rd_rinex_head.o) \
	$(TRACK_LIB)(rd_rinex_obs.o) \
	$(TRACK_LIB)(rx_obs.o) \
	$(TRACK_LIB)(openfile.o) \
	$(TRACK_LIB)(remap_rx.o) \
	$(TRACK_LIB)(read_nav.o) \
	$(TRACK_LIB)(smooth.o) \
	$(TRACK_LIB)(read_antmod.o) \
	$(TRACK_LIB)(intantmod.o) \
	$(TRACK_LIB)(ionex_lib.o) \
	$(TRACK_LIB)(cslip_repair.o) \
	$(TRACK_LIB)(assign_dattyp.o) \
	$(TRACK_LIB)(check_gnss.o) \
	$(TRACK_LIB)(NLx_est.o) \
	$(TRACK_LIB)(zero_vec8.o) \
	$(TRACK_LIB)(step2diu.o) \
	$(TRACK_LIB)(step2lon.o) \
	$(TRACK_LIB)(norm8.o) \
	$(TRACK_LIB)(sprod.o) \
	$(TRACK_LIB)(st1idiu.o) \
	$(TRACK_LIB)(st1isem.o) \
	$(TRACK_LIB)(st1l1.o) \
	$(TRACK_LIB)(dehanttideinel.o) \
	$(TRACK_LIB)(taiutc.o) \
	$(TRACK_LIB)(ephred.o) \
	$(TRACK_LIB)(ephtrp.o) \
	$(TRACK_LIB)(admint.o) \
	$(TRACK_LIB)(etutc.o) \
	$(TRACK_LIB)(eval.o) \
	$(TRACK_LIB)(juldat.o) \
	$(TRACK_LIB)(leap.o) \
	$(TRACK_LIB)(mday.o) \
	$(TRACK_LIB)(recurs.o) \
	$(TRACK_LIB)(shells.o) \
	$(TRACK_LIB)(spline.o) \
	$(TRACK_LIB)(toymd.o) \
	$(TRACK_LIB)(tdfrph.o) \
	$(TRACK_LIB)(hardisp.o) \
	$(TRACK_LIB)(vmf1_grid.o)
	ranlib track_lib.a

$(TRACK_LIB)(assign_con.o): makexk_cmds_bd.h
$(TRACK_LIB)(comp_clock.o): track_com.h
$(TRACK_LIB)(earth_to_inert.o): ../includes/const_param.h
$(TRACK_LIB)(get_elev.o): ../includes/const_param.h
$(TRACK_LIB)(get_miss.o): ../includes/const_param.h ../includes/xfile_def.h makexk_cmds_bd.h track_com.h
$(TRACK_LIB)(get_obsfile.o): makexk_cmds_bd.h track_com.h
$(TRACK_LIB)(get_site.o): ../includes/xfile_def.h makexk_cmds_bd.h track_com.h
$(TRACK_LIB)(get_svsrun.o): track_com.h
$(TRACK_LIB)(neph_to_xyz.o): track_com.h
$(TRACK_LIB)(read_bat.o): ../includes/xfile_def.h makexk_cmds_bd.h track_com.h
$(TRACK_LIB)(read_nav.o): track_com.h
$(TRACK_LIB)(read_nav_clk.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(save_head.o): ../includes/xfile_def.h track_com.h
$(TRACK_LIB)(clean_sp3_clk.o): track_com.h
$(TRACK_LIB)(init_track.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(read_all_rinex.o): ../includes/const_param.h ../includes/xfile_def.h track_com.h
$(TRACK_LIB)(est_pos.o): ../includes/const_param.h ../includes/xfile_def.h track_com.h vmf_com.h
$(TRACK_LIB)(theory.o): ../includes/const_param.h track_com.h vmf_com.h
$(TRACK_LIB)(search.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(atm_lib.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(rd_rinex_head.o): ../includes/xfile_def.h
$(TRACK_LIB)(rd_rinex_obs.o): ../includes/const_param.h ../includes/xfile_def.h track_com.h
$(TRACK_LIB)(rx_obs.o): ../includes/const_param.h ../includes/xfile_def.h
$(TRACK_LIB)(remap_rx.o): ../../libraries/includes/freq_def.h ../includes/xfile_def.h track_com.h
$(TRACK_LIB)(read_nav.o): track_com.h
$(TRACK_LIB)(smooth.o): track_com.h
$(TRACK_LIB)(read_antmod.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(ionex_lib.o): ../includes/const_param.h ../includes/xfile_def.h track_com.h
$(TRACK_LIB)(cslip_repair.o): ../includes/const_param.h track_com.h
$(TRACK_LIB)(assign_dattyp.o): ../includes/xfile_def.h
$(TRACK_LIB)(NLx_est.o): track_com.h
$(TRACK_LIB)(taiutc.o): track_com.h
$(TRACK_LIB)(ephred.o): track_com.h
$(TRACK_LIB)(hardisp.o): track_com.h
$(TRACK_LIB)(vmf1_grid.o): vmf_com.h

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
