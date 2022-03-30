# Makefile for utils on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:37:19 UTC 2021

all: apr_to_vel block corcom velrot merge_apr_vel merge_apr disp2apr doy ndoy edit_ext edit_apr edit_val gapr_to_l genrelref gen_sng gen_uselist gen_ssu getrel glbtog mit_to_iers plate cvframe ppb poles seasonal swaph unify_apr unify_svs velcom global_sel tgrep netsel phase_ext cortran stinf_to_rename memtest nutIAU2000 make_rms glist2cmd vel2stab sortrms grw lrw chi_to_rename fixstinf break_to_rename tdefnode2vel merge_renames



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

GEN_SSU_LIB = gen_ssu_lib.a
SLIB = ../Khandlers/KalObs_hndl_lib.a ../glinit/glinit_lib.a ../globk/globk_lib.a ../glfor/glfor_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../gamit/lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

apr_to_vel : apr_to_vel.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) apr_to_vel.f  $(SLIB) -o apr_to_vel
	rm -f apr_to_vel.o

block : block.f  $(SLIB) block.h
	$(FTN) $(FFLAGS) block.f  $(SLIB) -o block
	rm -f block.o

corcom : corcom.f  $(SLIB) ../includes/const_param.h corcom.h
	$(FTN) $(FFLAGS) corcom.f  $(SLIB) -o corcom
	rm -f corcom.o

velrot : velrot.f  $(SLIB) ../includes/const_param.h velrot.h
	$(FTN) $(FFLAGS) velrot.f  $(SLIB) -o velrot
	rm -f velrot.o

merge_apr_vel : merge_apr_vel.f  $(SLIB) ../includes/const_param.h merge_apr_vel.h
	$(FTN) $(FFLAGS) merge_apr_vel.f  $(SLIB) -o merge_apr_vel
	rm -f merge_apr_vel.o

merge_apr : merge_apr.f  $(SLIB) ../gen_util/nonsec_util.h
	$(FTN) $(FFLAGS) merge_apr.f  $(SLIB) -o merge_apr
	rm -f merge_apr.o

disp2apr : disp2apr.f  $(SLIB) 
	$(FTN) $(FFLAGS) disp2apr.f  $(SLIB) -o disp2apr
	rm -f disp2apr.o

doy : doy.f  $(SLIB) 
	$(FTN) $(FFLAGS) doy.f  $(SLIB) -o doy
	rm -f doy.o

ndoy : ndoy.f  $(SLIB) 
	$(FTN) $(FFLAGS) ndoy.f  $(SLIB) -o ndoy
	rm -f ndoy.o

edit_ext : edit_ext.f  $(SLIB) 
	$(FTN) $(FFLAGS) edit_ext.f  $(SLIB) -o edit_ext
	rm -f edit_ext.o

edit_apr : edit_apr.f  $(SLIB) 
	$(FTN) $(FFLAGS) edit_apr.f  $(SLIB) -o edit_apr
	rm -f edit_apr.o

edit_val : edit_val.f  $(SLIB) 
	$(FTN) $(FFLAGS) edit_val.f  $(SLIB) -o edit_val
	rm -f edit_val.o

gapr_to_l : gapr_to_l.f  $(SLIB) ../globk/globk_cmd_bd.f ../includes/const_param.h ../includes/globk_cmds.h ../includes/globk_common.h ../includes/kalman_param.h gapr_to_l.h
	$(FTN) $(FFLAGS) gapr_to_l.f  $(SLIB) -o gapr_to_l
	rm -f gapr_to_l.o

genrelref : genrelref.f  $(SLIB) 
	$(FTN) $(FFLAGS) genrelref.f  $(SLIB) -o genrelref
	rm -f genrelref.o

gen_sng : gen_sng.f  $(SLIB) 
	$(FTN) $(FFLAGS) gen_sng.f  $(SLIB) -o gen_sng
	rm -f gen_sng.o

gen_uselist : gen_uselist.f  $(SLIB) 
	$(FTN) $(FFLAGS) gen_uselist.f  $(SLIB) -o gen_uselist
	rm -f gen_uselist.o

gen_ssu : gen_ssu.f $(GEN_SSU_LIB) $(SLIB) ../includes/kalman_param.h ../includes/obs_header.h ../includes/sd_common.h
	$(FTN) $(FFLAGS) gen_ssu.f $(GEN_SSU_LIB) $(SLIB) -o gen_ssu
	rm -f gen_ssu.o

getrel : getrel.f  $(SLIB) ../includes/kalman_param.h
	$(FTN) $(FFLAGS) getrel.f  $(SLIB) -o getrel
	rm -f getrel.o

glbtog : glbtog.f  $(SLIB) glbtog.h
	$(FTN) $(FFLAGS) glbtog.f  $(SLIB) -o glbtog
	rm -f glbtog.o

mit_to_iers : mit_to_iers.f  $(SLIB) 
	$(FTN) $(FFLAGS) mit_to_iers.f  $(SLIB) -o mit_to_iers
	rm -f mit_to_iers.o

plate : plate.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) plate.f  $(SLIB) -o plate
	rm -f plate.o

cvframe : cvframe.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) cvframe.f  $(SLIB) -o cvframe
	rm -f cvframe.o

ppb : ppb.f  $(SLIB) ppb.h
	$(FTN) $(FFLAGS) ppb.f  $(SLIB) -o ppb
	rm -f ppb.o

poles : poles.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) poles.f  $(SLIB) -o poles
	rm -f poles.o

seasonal : seasonal.f  $(SLIB) ../includes/const_param.h seasonal.h
	$(FTN) $(FFLAGS) seasonal.f  $(SLIB) -o seasonal
	rm -f seasonal.o

swaph : swaph.f  $(SLIB) ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h ../includes/svinf_def.h
	$(FTN) $(FFLAGS) swaph.f  $(SLIB) -o swaph
	rm -f swaph.o

unify_apr : unify_apr.f  $(SLIB) 
	$(FTN) $(FFLAGS) unify_apr.f  $(SLIB) -o unify_apr
	rm -f unify_apr.o

unify_svs : unify_svs.f  $(SLIB) 
	$(FTN) $(FFLAGS) unify_svs.f  $(SLIB) -o unify_svs
	rm -f unify_svs.o

velcom : velcom.f  $(SLIB) velcom.h
	$(FTN) $(FFLAGS) velcom.f  $(SLIB) -o velcom
	rm -f velcom.o

global_sel : global_sel.f  $(SLIB) 
	$(FTN) $(FFLAGS) global_sel.f  $(SLIB) -o global_sel
	rm -f global_sel.o

tgrep : tgrep.f  $(SLIB) 
	$(FTN) $(FFLAGS) tgrep.f  $(SLIB) -o tgrep
	rm -f tgrep.o

netsel : netsel.f  $(SLIB) ../includes/const_param.h netsel.h
	$(FTN) $(FFLAGS) netsel.f  $(SLIB) -o netsel
	rm -f netsel.o

phase_ext : phase_ext.f  $(SLIB) 
	$(FTN) $(FFLAGS) phase_ext.f  $(SLIB) -o phase_ext
	rm -f phase_ext.o

cortran : cortran.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) cortran.f  $(SLIB) -o cortran
	rm -f cortran.o

stinf_to_rename : stinf_to_rename.f  $(SLIB) 
	$(FTN) $(FFLAGS) stinf_to_rename.f  $(SLIB) -o stinf_to_rename
	rm -f stinf_to_rename.o

memtest : memtest.f  $(SLIB) 
	$(FTN) $(FFLAGS) memtest.f  $(SLIB) -o memtest
	rm -f memtest.o

nutIAU2000 : nutIAU2000.f  $(SLIB) 
	$(FTN) $(FFLAGS) nutIAU2000.f  $(SLIB) -o nutIAU2000
	rm -f nutIAU2000.o

make_rms : make_rms.f  $(SLIB) 
	$(FTN) $(FFLAGS) make_rms.f  $(SLIB) -o make_rms
	rm -f make_rms.o

glist2cmd : glist2cmd.f  $(SLIB) glist2cmd.h
	$(FTN) $(FFLAGS) glist2cmd.f  $(SLIB) -o glist2cmd
	rm -f glist2cmd.o

vel2stab : vel2stab.f  $(SLIB) 
	$(FTN) $(FFLAGS) vel2stab.f  $(SLIB) -o vel2stab
	rm -f vel2stab.o

sortrms : sortrms.f  $(SLIB) 
	$(FTN) $(FFLAGS) sortrms.f  $(SLIB) -o sortrms
	rm -f sortrms.o

grw : grw.f  $(SLIB) 
	$(FTN) $(FFLAGS) grw.f  $(SLIB) -o grw
	rm -f grw.o

lrw : lrw.f  $(SLIB) 
	$(FTN) $(FFLAGS) lrw.f  $(SLIB) -o lrw
	rm -f lrw.o

chi_to_rename : chi_to_rename.f  $(SLIB) 
	$(FTN) $(FFLAGS) chi_to_rename.f  $(SLIB) -o chi_to_rename
	rm -f chi_to_rename.o

fixstinf : fixstinf.f  $(SLIB) 
	$(FTN) $(FFLAGS) fixstinf.f  $(SLIB) -o fixstinf
	rm -f fixstinf.o

break_to_rename : break_to_rename.f  $(SLIB) 
	$(FTN) $(FFLAGS) break_to_rename.f  $(SLIB) -o break_to_rename
	rm -f break_to_rename.o

tdefnode2vel : tdefnode2vel.f  $(SLIB) 
	$(FTN) $(FFLAGS) tdefnode2vel.f  $(SLIB) -o tdefnode2vel
	rm -f tdefnode2vel.o

merge_renames : merge_renames.f  $(SLIB) 
	$(FTN) $(FFLAGS) merge_renames.f  $(SLIB) -o merge_renames
	rm -f merge_renames.o

$(GEN_SSU_LIB) : \
	$(GEN_SSU_LIB)(sd_comp.o)
	ranlib gen_ssu_lib.a


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
