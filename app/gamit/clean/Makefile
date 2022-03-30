# Makefile for clean on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:17:09 UTC 2021

all: cview scandd showrms



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

CVIEW_LIB = cview_lib.a
SCAN_LIB = scan_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

cview : cview.f $(CVIEW_LIB) $(SLIB) ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
	$(FTN) $(FFLAGS) cview.f $(CVIEW_LIB) $(SLIB) -L$(X11LIB) -lX11 -o cview
	rm -f cview.o

scandd : scandd.f $(SCAN_LIB) $(SLIB) ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
	$(FTN) $(FFLAGS) scandd.f $(SCAN_LIB) $(SLIB) -L$(X11LIB) -lX11 -o scandd
	rm -f scandd.o

showrms : showrms.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) showrms.f  $(SLIB) -o showrms
	rm -f showrms.o

$(CVIEW_LIB) : \
	$(CVIEW_LIB)(addadj.o) \
	$(CVIEW_LIB)(allan.o) \
	$(CVIEW_LIB)(bound.o) \
	$(CVIEW_LIB)(atm_grad_part.o) \
	$(CVIEW_LIB)(check_wl.o) \
	$(CVIEW_LIB)(combo.o) \
	$(CVIEW_LIB)(cversn.o) \
	$(CVIEW_LIB)(clkdif.o) \
	$(CVIEW_LIB)(editor.o) \
	$(CVIEW_LIB)(ed_bias.o) \
	$(CVIEW_LIB)(drv_scn.o) \
	$(CVIEW_LIB)(dr_srch.o) \
	$(CVIEW_LIB)(dr_ptch.o) \
	$(CVIEW_LIB)(ed_data.o) \
	$(CVIEW_LIB)(ed_elim.o) \
	$(CVIEW_LIB)(ed_find.o) \
	$(CVIEW_LIB)(ed_line.o) \
	$(CVIEW_LIB)(ed_list.o) \
	$(CVIEW_LIB)(ed_mlbs.o) \
	$(CVIEW_LIB)(ed_mnu0.o) \
	$(CVIEW_LIB)(ed_mnu1.o) \
	$(CVIEW_LIB)(ed_mnu2.o) \
	$(CVIEW_LIB)(ed_mnu3.o) \
	$(CVIEW_LIB)(ed_move.o) \
	$(CVIEW_LIB)(ed_ptch.o) \
	$(CVIEW_LIB)(ed_save.o) \
	$(CVIEW_LIB)(ed_slip.o) \
	$(CVIEW_LIB)(ed_subs1.o) \
	$(CVIEW_LIB)(ed_subs2.o) \
	$(CVIEW_LIB)(evlmem.o) \
	$(CVIEW_LIB)(finish.o) \
	$(CVIEW_LIB)(getser.o) \
	$(CVIEW_LIB)(getfil.o) \
	$(CVIEW_LIB)(gj.o) \
	$(CVIEW_LIB)(gmsg_dumb.o) \
	$(CVIEW_LIB)(memcof.o) \
	$(CVIEW_LIB)(multip.o) \
	$(CVIEW_LIB)(nicer.o) \
	$(CVIEW_LIB)(nextob.o) \
	$(CVIEW_LIB)(patch.o) \
	$(CVIEW_LIB)(pointers.o) \
	$(CVIEW_LIB)(polyft.o) \
	$(CVIEW_LIB)(readc.o) \
	$(CVIEW_LIB)(readm.o) \
	$(CVIEW_LIB)(readx.o) \
	$(CVIEW_LIB)(rinex.o) \
	$(CVIEW_LIB)(set_freqs.o) \
	$(CVIEW_LIB)(smooth.o) \
	$(CVIEW_LIB)(vcopy.o) \
	$(CVIEW_LIB)(writec.o) \
	$(CVIEW_LIB)(write_vcview.o) \
	$(CVIEW_LIB)(get_color_ro.o) \
	$(CVIEW_LIB)(csubs_gftn.o) \
	$(CVIEW_LIB)(gsubs_x_gftn.o)
	ranlib cview_lib.a

$(SCAN_LIB) : \
	$(SCAN_LIB)(addadj.o) \
	$(SCAN_LIB)(atm_grad_part.o) \
	$(SCAN_LIB)(bound.o) \
	$(SCAN_LIB)(check_wl.o) \
	$(SCAN_LIB)(combo.o) \
	$(SCAN_LIB)(clkdif.o) \
	$(SCAN_LIB)(finish.o) \
	$(SCAN_LIB)(getfil.o) \
	$(SCAN_LIB)(drv_scn.o) \
	$(SCAN_LIB)(dr_srch.o) \
	$(SCAN_LIB)(dr_ptch.o) \
	$(SCAN_LIB)(getser.o) \
	$(SCAN_LIB)(gmsg_dumb.o) \
	$(SCAN_LIB)(multip.o) \
	$(SCAN_LIB)(patch.o) \
	$(SCAN_LIB)(nextob.o) \
	$(SCAN_LIB)(pointers.o) \
	$(SCAN_LIB)(readc.o) \
	$(SCAN_LIB)(readm.o) \
	$(SCAN_LIB)(readx.o) \
	$(SCAN_LIB)(rinex.o) \
	$(SCAN_LIB)(scanner.o) \
	$(SCAN_LIB)(vcopy.o) \
	$(SCAN_LIB)(writec.o) \
	$(SCAN_LIB)(write_vcview.o) \
	$(SCAN_LIB)(cversn.o) \
	$(SCAN_LIB)(get_color_ro.o) \
	$(SCAN_LIB)(csubs_gftn.o) \
	$(SCAN_LIB)(gsubs_x_gftn.o)
	ranlib scan_lib.a

$(CVIEW_LIB)(addadj.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(bound.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(check_wl.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(combo.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(editor.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_bias.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(drv_scn.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(dr_srch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(dr_ptch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(ed_data.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_elim.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_find.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_line.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_list.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_mlbs.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_mnu0.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_mnu1.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_mnu3.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_move.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_ptch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_save.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_slip.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(CVIEW_LIB)(ed_subs1.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(ed_subs2.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(finish.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(getser.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(getfil.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(nextob.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(patch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(pointers.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(readc.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(readm.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(readx.o): ../../libraries/includes/freq_def.h ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(rinex.o): ../../libraries/includes/freq_def.h ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(set_freqs.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(vcopy.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(CVIEW_LIB)(writec.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(CVIEW_LIB)(write_vcview.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(addadj.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(SCAN_LIB)(bound.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(check_wl.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(combo.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(finish.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(getfil.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(drv_scn.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(dr_srch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/macdep.h ../includes/makex.h
$(SCAN_LIB)(dr_ptch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(getser.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(patch.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(nextob.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(pointers.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(SCAN_LIB)(readc.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(readm.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(SCAN_LIB)(readx.o): ../../libraries/includes/freq_def.h ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(SCAN_LIB)(rinex.o): ../../libraries/includes/freq_def.h ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(scanner.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(vcopy.o): ../includes/cview.h ../includes/dimpar.h ../includes/makex.h
$(SCAN_LIB)(writec.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(SCAN_LIB)(write_vcview.o): ../includes/cview.h ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h

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
