# Makefile for utils on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:26:42 UTC 2021

all: addvec argo2rx argo2fic argo2nav atx2svnav autecl bsx2dcbtab check_lfile check_siteid chicurve cml confpol convert_antpcv conveuler convert_jplybias countx crx2rnx rnx2crx cvedt dcbtab dcbtab2 decyr_atm emerge emntab erptab eulervel fic2nav fixsst fixssi fixx gmt2vel interp_log_tab make_stnfolist merge_rinex metutil ne2azel poly01 poly01r plotk profile pull_bias rx2apr scan_rinex stnfo_cont svdlin svdsig vector_stat vel2hist vexclude wbslfilt testephred testnut



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

CONVEULER_LIB = conveuler_lib.a
EULERVEL_LIB = eulervel_lib.a
METUTIL_LIB = metutil_lib.a
VEXCLUDE_LIB = vexclude_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a ../../libraries/matrix/kinv_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

addvec : addvec.f  $(SLIB) 
	$(FTN) $(FFLAGS) addvec.f  $(SLIB) -o addvec
	rm -f addvec.o

argo2rx : argo2rx.f  $(SLIB) ../includes/dimpar.h ../includes/makex.h
	$(FTN) $(FFLAGS) argo2rx.f  $(SLIB) -o argo2rx
	rm -f argo2rx.o

argo2fic : argo2fic.f  $(SLIB) ../includes/argo2fic.h
	$(FTN) $(FFLAGS) argo2fic.f  $(SLIB) -o argo2fic
	rm -f argo2fic.o

argo2nav : argo2nav.f  $(SLIB) ../includes/argo2fic.h
	$(FTN) $(FFLAGS) argo2nav.f  $(SLIB) -o argo2nav
	rm -f argo2nav.o

atx2svnav : atx2svnav.f  $(SLIB) 
	$(FTN) $(FFLAGS) atx2svnav.f  $(SLIB) -o atx2svnav
	rm -f atx2svnav.o

autecl : autecl.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) autecl.f  $(SLIB) -o autecl
	rm -f autecl.o

bsx2dcbtab : bsx2dcbtab.f  $(SLIB) 
	$(FTN) $(FFLAGS) bsx2dcbtab.f  $(SLIB) -o bsx2dcbtab
	rm -f bsx2dcbtab.o

check_lfile : check_lfile.f  $(SLIB) 
	$(FTN) $(FFLAGS) check_lfile.f  $(SLIB) -o check_lfile
	rm -f check_lfile.o

check_siteid : check_siteid.f  $(SLIB) 
	$(FTN) $(FFLAGS) check_siteid.f  $(SLIB) -o check_siteid
	rm -f check_siteid.o

chicurve : chicurve.f  $(SLIB) 
	$(FTN) $(FFLAGS) chicurve.f  $(SLIB) -o chicurve
	rm -f chicurve.o

cml : cml.f  $(SLIB) 
	$(FTN) $(FFLAGS) cml.f  $(SLIB) -o cml
	rm -f cml.o

confpol : confpol.f  $(SLIB) 
	$(FTN) $(FFLAGS) confpol.f  $(SLIB) -o confpol
	rm -f confpol.o

convert_antpcv : convert_antpcv.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) convert_antpcv.f  $(SLIB) -o convert_antpcv
	rm -f convert_antpcv.o

conveuler : conveuler.f $(CONVEULER_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) conveuler.f $(CONVEULER_LIB) $(SLIB) -o conveuler
	rm -f conveuler.o

convert_jplybias : convert_jplybias.f  $(SLIB) 
	$(FTN) $(FFLAGS) convert_jplybias.f  $(SLIB) -o convert_jplybias
	rm -f convert_jplybias.o

countx : countx.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) countx.f  $(SLIB) -o countx
	rm -f countx.o

crx2rnx : crx2rnx.c  $(SLIB) 
	$(CC) $(CFLAGS) crx2rnx.c $(SLIB) -lm  -o crx2rnx
	rm -f crx2rnx.o

rnx2crx : rnx2crx.c  $(SLIB) 
	$(CC) $(CFLAGS) rnx2crx.c $(SLIB) -lm  -o rnx2crx
	rm -f rnx2crx.o

cvedt : cvedt.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) cvedt.f  $(SLIB) -o cvedt
	rm -f cvedt.o

dcbtab : dcbtab.f  $(SLIB) 
	$(FTN) $(FFLAGS) dcbtab.f  $(SLIB) -o dcbtab
	rm -f dcbtab.o

dcbtab2 : dcbtab2.f  $(SLIB) 
	$(FTN) $(FFLAGS) dcbtab2.f  $(SLIB) -o dcbtab2
	rm -f dcbtab2.o

decyr_atm : decyr_atm.f  $(SLIB) 
	$(FTN) $(FFLAGS) decyr_atm.f  $(SLIB) -o decyr_atm
	rm -f decyr_atm.o

emerge : emerge.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) emerge.f  $(SLIB) -o emerge
	rm -f emerge.o

emntab : emntab.f  $(SLIB) 
	$(FTN) $(FFLAGS) emntab.f  $(SLIB) -o emntab
	rm -f emntab.o

erptab : erptab.f  $(SLIB) 
	$(FTN) $(FFLAGS) erptab.f  $(SLIB) -o erptab
	rm -f erptab.o

eulervel : eulervel.f $(EULERVEL_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) eulervel.f $(EULERVEL_LIB) $(SLIB) -o eulervel
	rm -f eulervel.o

fic2nav : fic2nav.f  $(SLIB) ../includes/argo2fic.h ../includes/makex.h
	$(FTN) $(FFLAGS) fic2nav.f  $(SLIB) -o fic2nav
	rm -f fic2nav.o

fixsst : fixsst.f  $(SLIB) 
	$(FTN) $(FFLAGS) fixsst.f  $(SLIB) -o fixsst
	rm -f fixsst.o

fixssi : fixssi.f  $(SLIB) 
	$(FTN) $(FFLAGS) fixssi.f  $(SLIB) -o fixssi
	rm -f fixssi.o

fixx : fixx.f  $(SLIB) 
	$(FTN) $(FFLAGS) fixx.f  $(SLIB) -o fixx
	rm -f fixx.o

gmt2vel : gmt2vel.f  $(SLIB) 
	$(FTN) $(FFLAGS) gmt2vel.f  $(SLIB) -o gmt2vel
	rm -f gmt2vel.o

interp_log_tab : interp_log_tab.f  $(SLIB) 
	$(FTN) $(FFLAGS) interp_log_tab.f  $(SLIB) -o interp_log_tab
	rm -f interp_log_tab.o

make_stnfolist : make_stnfolist.f  $(SLIB) 
	$(FTN) $(FFLAGS) make_stnfolist.f  $(SLIB) -o make_stnfolist
	rm -f make_stnfolist.o

merge_rinex : merge_rinex.f  $(SLIB) ../includes/dimpar.h ../includes/makex.h
	$(FTN) $(FFLAGS) merge_rinex.f  $(SLIB) -o merge_rinex
	rm -f merge_rinex.o

metutil : metutil.f $(METUTIL_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) metutil.f $(METUTIL_LIB) $(SLIB) -o metutil
	rm -f metutil.o

ne2azel : ne2azel.f  $(SLIB) 
	$(FTN) $(FFLAGS) ne2azel.f  $(SLIB) -o ne2azel
	rm -f ne2azel.o

poly01 : poly01.f  $(SLIB) 
	$(FTN) $(FFLAGS) poly01.f  $(SLIB) -o poly01
	rm -f poly01.o

poly01r : poly01r.f  $(SLIB) 
	$(FTN) $(FFLAGS) poly01r.f  $(SLIB) -o poly01r
	rm -f poly01r.o

plotk : plotk.f  $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) plotk.f  $(SLIB) -o plotk
	rm -f plotk.o

profile : profile.f  $(SLIB) 
	$(FTN) $(FFLAGS) profile.f  $(SLIB) -o profile
	rm -f profile.o

pull_bias : pull_bias.f  $(SLIB) 
	$(FTN) $(FFLAGS) pull_bias.f  $(SLIB) -o pull_bias
	rm -f pull_bias.o

rx2apr : rx2apr.f  $(SLIB) 
	$(FTN) $(FFLAGS) rx2apr.f  $(SLIB) -o rx2apr
	rm -f rx2apr.o

scan_rinex : scan_rinex.f  $(SLIB) ../includes/dimpar.h rinex.h
	$(FTN) $(FFLAGS) scan_rinex.f  $(SLIB) -o scan_rinex
	rm -f scan_rinex.o

stnfo_cont : stnfo_cont.f  $(SLIB) 
	$(FTN) $(FFLAGS) stnfo_cont.f  $(SLIB) -o stnfo_cont
	rm -f stnfo_cont.o

svdlin : svdlin.f  $(SLIB) 
	$(FTN) $(FFLAGS) svdlin.f  $(SLIB) -o svdlin
	rm -f svdlin.o

svdsig : svdsig.f  $(SLIB) 
	$(FTN) $(FFLAGS) svdsig.f  $(SLIB) -o svdsig
	rm -f svdsig.o

vector_stat : vector_stat.f  $(SLIB) 
	$(FTN) $(FFLAGS) vector_stat.f  $(SLIB) -o vector_stat
	rm -f vector_stat.o

vel2hist : vel2hist.f  $(SLIB) 
	$(FTN) $(FFLAGS) vel2hist.f  $(SLIB) -o vel2hist
	rm -f vel2hist.o

vexclude : vexclude.f $(VEXCLUDE_LIB) $(SLIB) ../includes/dimpar.h
	$(FTN) $(FFLAGS) vexclude.f $(VEXCLUDE_LIB) $(SLIB) -o vexclude
	rm -f vexclude.o

wbslfilt : wbslfilt.f  $(SLIB) 
	$(FTN) $(FFLAGS) wbslfilt.f  $(SLIB) -o wbslfilt
	rm -f wbslfilt.o

testephred : testephred.f  $(SLIB) ../includes/arc.h ../includes/dimpar.h ../includes/global.h ../includes/units.h
	$(FTN) $(FFLAGS) testephred.f  $(SLIB) -o testephred
	rm -f testephred.o

testnut : testnut.f  $(SLIB) 
	$(FTN) $(FFLAGS) testnut.f  $(SLIB) -o testnut
	rm -f testnut.o

$(CONVEULER_LIB) : \
	$(CONVEULER_LIB)(xyzpole.o) \
	$(CONVEULER_LIB)(velpole.o)
	ranlib conveuler_lib.a

$(EULERVEL_LIB) : \
	$(EULERVEL_LIB)(xyzpole.o) \
	$(EULERVEL_LIB)(velpole.o)
	ranlib eulervel_lib.a

$(METUTIL_LIB) : \
	$(METUTIL_LIB)(uversn.o) \
	$(METUTIL_LIB)(geoc_to_geod.o)
	ranlib metutil_lib.a

$(VEXCLUDE_LIB) : \
	$(VEXCLUDE_LIB)(read_q.o) \
	$(VEXCLUDE_LIB)(change_excl.o) \
	$(VEXCLUDE_LIB)(uversn.o)
	ranlib vexclude_lib.a

$(VEXCLUDE_LIB)(read_q.o): ../includes/dimpar.h
$(VEXCLUDE_LIB)(change_excl.o): ../includes/dimpar.h

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
