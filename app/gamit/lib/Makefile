# Makefile for lib on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:13:11 UTC 2021

all: gamit_lib



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

GAMIT_LIB = gamit_lib.a
SLIB = ../../libraries/comlib/com_lib.a ./libsofa.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

gamit_lib : $(GAMIT_LIB)

$(GAMIT_LIB) : \
	$(GAMIT_LIB)(amag3.o) \
	$(GAMIT_LIB)(angvec.o) \
	$(GAMIT_LIB)(ant_alias.o) \
	$(GAMIT_LIB)(argmnt.o) \
	$(GAMIT_LIB)(assign_srpnames.o) \
	$(GAMIT_LIB)(blank.o) \
	$(GAMIT_LIB)(block_name.o) \
	$(GAMIT_LIB)(cdate.o) \
	$(GAMIT_LIB)(check.o) \
	$(GAMIT_LIB)(check_y2k.o) \
	$(GAMIT_LIB)(clkera.o) \
	$(GAMIT_LIB)(copens.o) \
	$(GAMIT_LIB)(count_arg.o) \
	$(GAMIT_LIB)(crd_file_type.o) \
	$(GAMIT_LIB)(cross.o) \
	$(GAMIT_LIB)(cross_unit.o) \
	$(GAMIT_LIB)(crs_trs.o) \
	$(GAMIT_LIB)(dayjul.o) \
	$(GAMIT_LIB)(daynum.o) \
	$(GAMIT_LIB)(dblas1.o) \
	$(GAMIT_LIB)(dblas2.o) \
	$(GAMIT_LIB)(dblas3.o) \
	$(GAMIT_LIB)(dlapack.o) \
	$(GAMIT_LIB)(dlinpk.o) \
	$(GAMIT_LIB)(decyrs.o) \
	$(GAMIT_LIB)(degdms.o) \
	$(GAMIT_LIB)(dmsdeg.o) \
	$(GAMIT_LIB)(dmsrad.o) \
	$(GAMIT_LIB)(dot.o) \
	$(GAMIT_LIB)(doygwk.o) \
	$(GAMIT_LIB)(dvdot.o) \
	$(GAMIT_LIB)(dvswp.o) \
	$(GAMIT_LIB)(ephdrd.o) \
	$(GAMIT_LIB)(errcod.o) \
	$(GAMIT_LIB)(even_minute.o) \
	$(GAMIT_LIB)(evrtcf.o) \
	$(GAMIT_LIB)(fcheck.o) \
	$(GAMIT_LIB)(fix_y2k.o) \
	$(GAMIT_LIB)(frac.o) \
	$(GAMIT_LIB)(funarg.o) \
	$(GAMIT_LIB)(funcof.o) \
	$(GAMIT_LIB)(geninv1.o) \
	$(GAMIT_LIB)(geoxyz.o) \
	$(GAMIT_LIB)(geodms_xyz.o) \
	$(GAMIT_LIB)(getdat.o) \
	$(GAMIT_LIB)(getcmd.o) \
	$(GAMIT_LIB)(gettim.o) \
	$(GAMIT_LIB)(get_iut1pol.o) \
	$(GAMIT_LIB)(hfeop_util.o) \
	$(GAMIT_LIB)(hfeop_desai.o) \
	$(GAMIT_LIB)(hfeop_gipson.o) \
	$(GAMIT_LIB)(getusr.o) \
	$(GAMIT_LIB)(get_gnss.o) \
	$(GAMIT_LIB)(gsatel.o) \
	$(GAMIT_LIB)(gpsblock.o) \
	$(GAMIT_LIB)(gpt.o) \
	$(GAMIT_LIB)(gpt2.o) \
	$(GAMIT_LIB)(gpt3.o) \
	$(GAMIT_LIB)(hisub.o) \
	$(GAMIT_LIB)(iarray.o) \
	$(GAMIT_LIB)(ic4array.o) \
	$(GAMIT_LIB)(iclarg.o) \
	$(GAMIT_LIB)(idoy.o) \
	$(GAMIT_LIB)(ier2df.o) \
	$(GAMIT_LIB)(ierrcm.o) \
	$(GAMIT_LIB)(imenu.o) \
	$(GAMIT_LIB)(indexi.o) \
	$(GAMIT_LIB)(indexx.o) \
	$(GAMIT_LIB)(interp_azel.o) \
	$(GAMIT_LIB)(inver2.o) \
	$(GAMIT_LIB)(itimdif.o) \
	$(GAMIT_LIB)(jel.o) \
	$(GAMIT_LIB)(jelf.o) \
	$(GAMIT_LIB)(julday.o) \
	$(GAMIT_LIB)(jul2hms.o) \
	$(GAMIT_LIB)(lask.o) \
	$(GAMIT_LIB)(lbias.o) \
	$(GAMIT_LIB)(lgood.o) \
	$(GAMIT_LIB)(lift_arg.o) \
	$(GAMIT_LIB)(linear.o) \
	$(GAMIT_LIB)(lloel.o) \
	$(GAMIT_LIB)(lmarg.o) \
	$(GAMIT_LIB)(lower1.o) \
	$(GAMIT_LIB)(lowerc.o) \
	$(GAMIT_LIB)(lowers.o) \
	$(GAMIT_LIB)(lread.o) \
	$(GAMIT_LIB)(lversn.o) \
	$(GAMIT_LIB)(lunred.o) \
	$(GAMIT_LIB)(mapamp.o) \
	$(GAMIT_LIB)(match_name.o) \
	$(GAMIT_LIB)(matmpy.o) \
	$(GAMIT_LIB)(mchkey.o) \
	$(GAMIT_LIB)(MHB_2000.o) \
	$(GAMIT_LIB)(mmply.o) \
	$(GAMIT_LIB)(monday.o) \
	$(GAMIT_LIB)(mopens.o) \
	$(GAMIT_LIB)(nblen.o) \
	$(GAMIT_LIB)(newchr.o) \
	$(GAMIT_LIB)(normalise.o) \
	$(GAMIT_LIB)(nutred.o) \
	$(GAMIT_LIB)(nuttab.o) \
	$(GAMIT_LIB)(nydays.o) \
	$(GAMIT_LIB)(ocearg.o) \
	$(GAMIT_LIB)(ocearg2.o) \
	$(GAMIT_LIB)(otlcmc.o) \
	$(GAMIT_LIB)(pjdhms.o) \
	$(GAMIT_LIB)(pnrot.o) \
	$(GAMIT_LIB)(pns.o) \
	$(GAMIT_LIB)(polred.o) \
	$(GAMIT_LIB)(prces.o) \
	$(GAMIT_LIB)(printmat.o) \
	$(GAMIT_LIB)(proper.o) \
	$(GAMIT_LIB)(raddms.o) \
	$(GAMIT_LIB)(rand5.o) \
	$(GAMIT_LIB)(rdsest.o) \
	$(GAMIT_LIB)(read_antmod.o) \
	$(GAMIT_LIB)(read_antex.o) \
	$(GAMIT_LIB)(read_antex_head.o) \
	$(GAMIT_LIB)(readc1.o) \
	$(GAMIT_LIB)(readc2.o) \
	$(GAMIT_LIB)(readc3.o) \
	$(GAMIT_LIB)(readc4.o) \
	$(GAMIT_LIB)(readc5.o) \
	$(GAMIT_LIB)(readdf.o) \
	$(GAMIT_LIB)(reade.o) \
	$(GAMIT_LIB)(readj.o) \
	$(GAMIT_LIB)(readm1.o) \
	$(GAMIT_LIB)(readm2.o) \
	$(GAMIT_LIB)(readm3.o) \
	$(GAMIT_LIB)(read_gdatum.o) \
	$(GAMIT_LIB)(read_gfile.o) \
	$(GAMIT_LIB)(read_rcvant.o) \
	$(GAMIT_LIB)(rfica.o) \
	$(GAMIT_LIB)(rixhed.o) \
	$(GAMIT_LIB)(rjust.o) \
	$(GAMIT_LIB)(rnhalf.o) \
	$(GAMIT_LIB)(rnqrtr.o) \
	$(GAMIT_LIB)(rotate_geod.o) \
	$(GAMIT_LIB)(rot_gfile.o) \
	$(GAMIT_LIB)(rotmat.o) \
	$(GAMIT_LIB)(rotsnp.o) \
	$(GAMIT_LIB)(rotsnp_sofa.o) \
	$(GAMIT_LIB)(round6h.o) \
	$(GAMIT_LIB)(rrinex.o) \
	$(GAMIT_LIB)(rrxhed.o) \
	$(GAMIT_LIB)(rsesfo.o) \
	$(GAMIT_LIB)(rsp3hd.o) \
	$(GAMIT_LIB)(rstnfo.o) \
	$(GAMIT_LIB)(runtim.o) \
	$(GAMIT_LIB)(satatt.o) \
	$(GAMIT_LIB)(secdif.o) \
	$(GAMIT_LIB)(secsum.o) \
	$(GAMIT_LIB)(sel_obtyp.o) \
	$(GAMIT_LIB)(shadow1.o) \
	$(GAMIT_LIB)(sidmat.o) \
	$(GAMIT_LIB)(sidtim.o) \
	$(GAMIT_LIB)(snp.o) \
	$(GAMIT_LIB)(solred.o) \
	$(GAMIT_LIB)(sort_string.o) \
	$(GAMIT_LIB)(sph2xyz.o) \
	$(GAMIT_LIB)(srotat.o) \
	$(GAMIT_LIB)(svnav_read.o) \
	$(GAMIT_LIB)(suicid.o) \
	$(GAMIT_LIB)(thdred.o) \
	$(GAMIT_LIB)(strjst.o) \
	$(GAMIT_LIB)(taiutc.o) \
	$(GAMIT_LIB)(timcon.o) \
	$(GAMIT_LIB)(timdif.o) \
	$(GAMIT_LIB)(time_TT.o) \
	$(GAMIT_LIB)(time_TT_sec.o) \
	$(GAMIT_LIB)(timinc.o) \
	$(GAMIT_LIB)(topens.o) \
	$(GAMIT_LIB)(transp.o) \
	$(GAMIT_LIB)(upnam1.o) \
	$(GAMIT_LIB)(upnam2.o) \
	$(GAMIT_LIB)(upnam3.o) \
	$(GAMIT_LIB)(upper1.o) \
	$(GAMIT_LIB)(upperc.o) \
	$(GAMIT_LIB)(uppers.o) \
	$(GAMIT_LIB)(utc2gps.o) \
	$(GAMIT_LIB)(ut1red.o) \
	$(GAMIT_LIB)(ut1tid.o) \
	$(GAMIT_LIB)(utlibr.o) \
	$(GAMIT_LIB)(wdms.o) \
	$(GAMIT_LIB)(wfica.o) \
	$(GAMIT_LIB)(writc1.o) \
	$(GAMIT_LIB)(writc2.o) \
	$(GAMIT_LIB)(writc3.o) \
	$(GAMIT_LIB)(writc4.o) \
	$(GAMIT_LIB)(writc5.o) \
	$(GAMIT_LIB)(writm1.o) \
	$(GAMIT_LIB)(writm2.o) \
	$(GAMIT_LIB)(writm3.o) \
	$(GAMIT_LIB)(wrxhed.o) \
	$(GAMIT_LIB)(wrxnav.o) \
	$(GAMIT_LIB)(wstnfo.o) \
	$(GAMIT_LIB)(wtime.o) \
	$(GAMIT_LIB)(wwrite.o) \
	$(GAMIT_LIB)(xyz2sph.o) \
	$(GAMIT_LIB)(xhdred.o) \
	$(GAMIT_LIB)(yaw_attit.o) \
	$(GAMIT_LIB)(ephred.o) \
	$(GAMIT_LIB)(ephtrp.o) \
	$(GAMIT_LIB)(read_svsinex.o)
	ranlib gamit_lib.a

$(GAMIT_LIB)(copens.o): ../includes/dimpar.h
$(GAMIT_LIB)(ephdrd.o): ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
$(GAMIT_LIB)(errcod.o): ../includes/errflg.h
$(GAMIT_LIB)(evrtcf.o): ../includes/arc.h ../includes/dimpar.h
$(GAMIT_LIB)(get_gnss.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(gsatel.o): ../includes/dimpar.h
$(GAMIT_LIB)(hisub.o): ../includes/dimpar.h
$(GAMIT_LIB)(ier2df.o): ../includes/errflg.h
$(GAMIT_LIB)(ierrcm.o): ../includes/errflg.h
$(GAMIT_LIB)(interp_azel.o): ../includes/dimpar.h
$(GAMIT_LIB)(inver2.o): ../includes/dimpar.h
$(GAMIT_LIB)(lbias.o): ../includes/errflg.h
$(GAMIT_LIB)(lgood.o): ../includes/errflg.h
$(GAMIT_LIB)(linear.o): ../includes/dimpar.h
$(GAMIT_LIB)(lloel.o): ../includes/errflg.h
$(GAMIT_LIB)(lmarg.o): ../includes/errflg.h
$(GAMIT_LIB)(lread.o): ../includes/dimpar.h ../includes/model.h ../includes/units.h
$(GAMIT_LIB)(lunred.o): ../includes/arc.h ../includes/dimpar.h ../includes/units.h
$(GAMIT_LIB)(mapamp.o): ../includes/makex.h
$(GAMIT_LIB)(MHB_2000.o): ../../kf/includes/const_param.h
$(GAMIT_LIB)(read_antmod.o): ../includes/dimpar.h
$(GAMIT_LIB)(read_antex.o): ../includes/dimpar.h
$(GAMIT_LIB)(read_antex_head.o): ../includes/dimpar.h
$(GAMIT_LIB)(readc1.o): ../includes/dimpar.h
$(GAMIT_LIB)(readc2.o): ../includes/dimpar.h
$(GAMIT_LIB)(readc3.o): ../includes/dimpar.h
$(GAMIT_LIB)(readc4.o): ../includes/dimpar.h
$(GAMIT_LIB)(readc5.o): ../includes/dimpar.h
$(GAMIT_LIB)(readdf.o): ../includes/dimpar.h
$(GAMIT_LIB)(reade.o): ../includes/makex.h
$(GAMIT_LIB)(readj.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(readm1.o): ../includes/dimpar.h
$(GAMIT_LIB)(readm2.o): ../includes/dimpar.h
$(GAMIT_LIB)(readm3.o): ../includes/dimpar.h
$(GAMIT_LIB)(read_gfile.o): ../includes/dimpar.h
$(GAMIT_LIB)(rixhed.o): ../includes/dimpar.h
$(GAMIT_LIB)(rot_gfile.o): ../includes/dimpar.h
$(GAMIT_LIB)(rrinex.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(rrxhed.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(rsesfo.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(rsp3hd.o): ../includes/dimpar.h
$(GAMIT_LIB)(satatt.o): ../includes/dimpar.h
$(GAMIT_LIB)(sel_obtyp.o): ../includes/makex.h
$(GAMIT_LIB)(shadow1.o): ../includes/dimpar.h
$(GAMIT_LIB)(solred.o): ../includes/arc.h ../includes/dimpar.h ../includes/units.h
$(GAMIT_LIB)(thdred.o): ../includes/dimpar.h
$(GAMIT_LIB)(writc1.o): ../includes/dimpar.h
$(GAMIT_LIB)(writc2.o): ../includes/dimpar.h
$(GAMIT_LIB)(writc3.o): ../includes/dimpar.h
$(GAMIT_LIB)(writc4.o): ../includes/dimpar.h
$(GAMIT_LIB)(writc5.o): ../../libraries/includes/const_param.h ../includes/dimpar.h
$(GAMIT_LIB)(writm1.o): ../includes/dimpar.h
$(GAMIT_LIB)(writm2.o): ../includes/dimpar.h
$(GAMIT_LIB)(writm3.o): ../includes/dimpar.h
$(GAMIT_LIB)(wrxhed.o): ../includes/dimpar.h ../includes/makex.h
$(GAMIT_LIB)(xhdred.o): ../includes/dimpar.h

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
