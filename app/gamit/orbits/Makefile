# Makefile for orbits on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:22:26 UTC 2021

all: bctot gtog orbfit orbdif ngstot tmerge ttoasc ttog ttoics ttot ttongs yawtab ytoasc y2orbex



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

BCTOT_LIB = bctot_lib.a
GTOG_LIB = gtog_lib.a
ORBFIT_LIB = orbfit_lib.a
ORBDIF_LIB = orbdif_lib.a
NGSTOT_LIB = ngstot_lib.a
TMERGE_LIB = tmerge_lib.a
TTOASC_LIB = ttoasc_lib.a
TTOG_LIB = ttog_lib.a
TTOICS_LIB = ttoics_lib.a
TTOT_LIB = ttot_lib.a
TTONGS_LIB = ttongs_lib.a
YAWTAB_LIB = yawtab_lib.a
YTOASC_LIB = ytoasc_lib.a
Y2ORBEX_LIB = y2orbex_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a ../lib/libsofa.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

bctot : bctot.f $(BCTOT_LIB) $(SLIB) ../includes/dimpar.h ../includes/orbits.h
	$(FTN) $(FFLAGS) bctot.f $(BCTOT_LIB) $(SLIB) -o bctot
	rm -f bctot.o

gtog : gtog.f $(GTOG_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) gtog.f $(GTOG_LIB) $(SLIB) -o gtog
	rm -f gtog.o

orbfit : orbfit.f $(ORBFIT_LIB) $(SLIB) ../includes/dimpar.h ../includes/orbits.h orbfit.h
	$(FTN) $(FFLAGS) orbfit.f $(ORBFIT_LIB) $(SLIB) -o orbfit
	rm -f orbfit.o

orbdif : orbdif.f $(ORBDIF_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) orbdif.f $(ORBDIF_LIB) $(SLIB) -o orbdif
	rm -f orbdif.o

ngstot : ngstot.f $(NGSTOT_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ngstot.f $(NGSTOT_LIB) $(SLIB) -o ngstot
	rm -f ngstot.o

tmerge : tmerge.f $(TMERGE_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) tmerge.f $(TMERGE_LIB) $(SLIB) -o tmerge
	rm -f tmerge.o

ttoasc : ttoasc.f $(TTOASC_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ttoasc.f $(TTOASC_LIB) $(SLIB) -o ttoasc
	rm -f ttoasc.o

ttog : ttog.f $(TTOG_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ttog.f $(TTOG_LIB) $(SLIB) -o ttog
	rm -f ttog.o

ttoics : ttoics.f $(TTOICS_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ttoics.f $(TTOICS_LIB) $(SLIB) -o ttoics
	rm -f ttoics.o

ttot : ttot.f $(TTOT_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ttot.f $(TTOT_LIB) $(SLIB) -o ttot
	rm -f ttot.o

ttongs : ttongs.f $(TTONGS_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ttongs.f $(TTONGS_LIB) $(SLIB) -o ttongs
	rm -f ttongs.o

yawtab : yawtab.f $(YAWTAB_LIB) $(SLIB) ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/orbits.h ../includes/units.h
	$(FTN) $(FFLAGS) yawtab.f $(YAWTAB_LIB) $(SLIB) -o yawtab
	rm -f yawtab.o

ytoasc : ytoasc.f $(YTOASC_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) ytoasc.f $(YTOASC_LIB) $(SLIB) -o ytoasc
	rm -f ytoasc.o

y2orbex : y2orbex.f $(Y2ORBEX_LIB) $(SLIB) ../includes/dimpar.h orbex.h
	$(FTN) $(FFLAGS) y2orbex.f $(Y2ORBEX_LIB) $(SLIB) -o y2orbex
	rm -f y2orbex.o

$(BCTOT_LIB) : \
	$(BCTOT_LIB)(arcinp.o) \
	$(BCTOT_LIB)(brdxyz.o) \
	$(BCTOT_LIB)(chksvs.o) \
	$(BCTOT_LIB)(closeb.o) \
	$(BCTOT_LIB)(closest_epoch.o) \
	$(BCTOT_LIB)(filnmb.o) \
	$(BCTOT_LIB)(gmake.o) \
	$(BCTOT_LIB)(openb.o) \
	$(BCTOT_LIB)(oversn.o) \
	$(BCTOT_LIB)(rotcrd.o) \
	$(BCTOT_LIB)(sortbc.o) \
	$(BCTOT_LIB)(spanup.o) \
	$(BCTOT_LIB)(thdrit.o) \
	$(BCTOT_LIB)(trot.o)
	ranlib bctot_lib.a

$(GTOG_LIB) : \
	$(GTOG_LIB)(open_eop.o) \
	$(GTOG_LIB)(oversn.o)
	ranlib gtog_lib.a

$(ORBFIT_LIB) : \
	$(ORBFIT_LIB)(elem.o) \
	$(ORBFIT_LIB)(eopart.o) \
	$(ORBFIT_LIB)(gpgt.o) \
	$(ORBFIT_LIB)(kepler.o) \
	$(ORBFIT_LIB)(norminc.o) \
	$(ORBFIT_LIB)(norm_solve.o) \
	$(ORBFIT_LIB)(oversn.o) \
	$(ORBFIT_LIB)(partl.o) \
	$(ORBFIT_LIB)(plt_postfit.o) \
	$(ORBFIT_LIB)(read_input.o) \
	$(ORBFIT_LIB)(tran_part.o) \
	$(ORBFIT_LIB)(write_g.o) \
	$(ORBFIT_LIB)(write_summary.o) \
	$(ORBFIT_LIB)(xyz2rac.o)
	ranlib orbfit_lib.a

$(ORBDIF_LIB) : \
	$(ORBDIF_LIB)(orbrms.o) \
	$(ORBDIF_LIB)(oversn.o) \
	$(ORBDIF_LIB)(ttemp.o)
	ranlib orbdif_lib.a

$(NGSTOT_LIB) : \
	$(NGSTOT_LIB)(closeb.o) \
	$(NGSTOT_LIB)(gmake.o) \
	$(NGSTOT_LIB)(getics.o) \
	$(NGSTOT_LIB)(ngsopn.o) \
	$(NGSTOT_LIB)(oversn.o) \
	$(NGSTOT_LIB)(openb.o) \
	$(NGSTOT_LIB)(rsp1hd.o) \
	$(NGSTOT_LIB)(rotcrd.o) \
	$(NGSTOT_LIB)(thdrit.o) \
	$(NGSTOT_LIB)(tdtrit.o) \
	$(NGSTOT_LIB)(trot.o)
	ranlib ngstot_lib.a

$(TMERGE_LIB) : \
	$(TMERGE_LIB)(thdrit.o) \
	$(TMERGE_LIB)(oversn.o)
	ranlib tmerge_lib.a

$(TTOASC_LIB) : \
	$(TTOASC_LIB)(oversn.o)
	ranlib ttoasc_lib.a

$(TTOG_LIB) : \
	$(TTOG_LIB)(gmake.o) \
	$(TTOG_LIB)(oversn.o)
	ranlib ttog_lib.a

$(TTOICS_LIB) : \
	$(TTOICS_LIB)(oversn.o)
	ranlib ttoics_lib.a

$(TTOT_LIB) : \
	$(TTOT_LIB)(closeb.o) \
	$(TTOT_LIB)(openb.o) \
	$(TTOT_LIB)(oversn.o) \
	$(TTOT_LIB)(rotcrd.o) \
	$(TTOT_LIB)(thdrit.o) \
	$(TTOT_LIB)(trot.o)
	ranlib ttot_lib.a

$(TTONGS_LIB) : \
	$(TTONGS_LIB)(closeb.o) \
	$(TTONGS_LIB)(openb.o) \
	$(TTONGS_LIB)(oversn.o) \
	$(TTONGS_LIB)(rdsvclk.o) \
	$(TTONGS_LIB)(rotcrd.o) \
	$(TTONGS_LIB)(sdtrit.o) \
	$(TTONGS_LIB)(trot.o) \
	$(TTONGS_LIB)(thdrit.o) \
	$(TTONGS_LIB)(wsp1hd.o) \
	$(TTONGS_LIB)(wsp3hd.o)
	ranlib ttongs_lib.a

$(YAWTAB_LIB) : \
	$(YAWTAB_LIB)(irnss_yaw.o) \
	$(YAWTAB_LIB)(kouba_yaw.o) \
	$(YAWTAB_LIB)(kouba_gps.o) \
	$(YAWTAB_LIB)(kouba_glonass.o) \
	$(YAWTAB_LIB)(kouba_beidou.o) \
	$(YAWTAB_LIB)(wang_beidou.o) \
	$(YAWTAB_LIB)(kouba_galileo.o) \
	$(YAWTAB_LIB)(mu_angle.o) \
	$(YAWTAB_LIB)(oversn.o) \
	$(YAWTAB_LIB)(svsun_angles.o)
	ranlib yawtab_lib.a

$(YTOASC_LIB) : \
	$(YTOASC_LIB)(oversn.o)
	ranlib ytoasc_lib.a

$(Y2ORBEX_LIB) : \
	$(Y2ORBEX_LIB)(oversn.o) \
	$(Y2ORBEX_LIB)(get_svcoords.o) \
	$(Y2ORBEX_LIB)(quaternion.o) \
	$(Y2ORBEX_LIB)(rotcrd.o) \
	$(Y2ORBEX_LIB)(write_orbex.o)
	ranlib y2orbex_lib.a

$(BCTOT_LIB)(arcinp.o): ../includes/dimpar.h
$(BCTOT_LIB)(chksvs.o): ../includes/dimpar.h
$(BCTOT_LIB)(closest_epoch.o): ../includes/dimpar.h ../includes/orbits.h
$(BCTOT_LIB)(gmake.o): ../includes/dimpar.h
$(BCTOT_LIB)(sortbc.o): ../includes/dimpar.h ../includes/orbits.h
$(BCTOT_LIB)(thdrit.o): ../includes/dimpar.h
$(BCTOT_LIB)(trot.o): ../includes/dimpar.h
$(ORBFIT_LIB)(eopart.o): ../includes/dimpar.h
$(ORBFIT_LIB)(norminc.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(norm_solve.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(partl.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(plt_postfit.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(read_input.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(tran_part.o): ../includes/dimpar.h ../includes/orbits.h
$(ORBFIT_LIB)(write_g.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBFIT_LIB)(write_summary.o): ../includes/dimpar.h ../includes/orbits.h orbfit.h
$(ORBDIF_LIB)(orbrms.o): ../includes/dimpar.h
$(ORBDIF_LIB)(ttemp.o): ../includes/dimpar.h
$(NGSTOT_LIB)(gmake.o): ../includes/dimpar.h
$(NGSTOT_LIB)(getics.o): ../includes/dimpar.h
$(NGSTOT_LIB)(rsp1hd.o): ../includes/dimpar.h
$(NGSTOT_LIB)(thdrit.o): ../includes/dimpar.h
$(NGSTOT_LIB)(tdtrit.o): ../includes/dimpar.h
$(NGSTOT_LIB)(trot.o): ../includes/dimpar.h
$(TMERGE_LIB)(thdrit.o): ../includes/dimpar.h
$(TTOG_LIB)(gmake.o): ../includes/dimpar.h
$(TTOT_LIB)(thdrit.o): ../includes/dimpar.h
$(TTOT_LIB)(trot.o): ../includes/dimpar.h
$(TTONGS_LIB)(rdsvclk.o): ../includes/dimpar.h
$(TTONGS_LIB)(sdtrit.o): ../includes/dimpar.h
$(TTONGS_LIB)(trot.o): ../includes/dimpar.h
$(TTONGS_LIB)(thdrit.o): ../includes/dimpar.h
$(TTONGS_LIB)(wsp1hd.o): ../includes/dimpar.h
$(TTONGS_LIB)(wsp3hd.o): ../includes/dimpar.h
$(YAWTAB_LIB)(irnss_yaw.o): ../includes/dimpar.h ../includes/orbits.h
$(YAWTAB_LIB)(kouba_yaw.o): ../includes/dimpar.h ../includes/orbits.h
$(YAWTAB_LIB)(kouba_gps.o): ../includes/dimpar.h ../includes/orbits.h
$(YAWTAB_LIB)(kouba_glonass.o): ../includes/dimpar.h ../includes/orbits.h
$(YAWTAB_LIB)(kouba_beidou.o): ../includes/dimpar.h ../includes/orbits.h
$(YAWTAB_LIB)(kouba_galileo.o): ../includes/dimpar.h
$(Y2ORBEX_LIB)(get_svcoords.o): ../includes/dimpar.h orbex.h
$(Y2ORBEX_LIB)(quaternion.o): ../includes/dimpar.h orbex.h
$(Y2ORBEX_LIB)(write_orbex.o): ../includes/dimpar.h orbex.h

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
