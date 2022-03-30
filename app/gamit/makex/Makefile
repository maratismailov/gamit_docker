# Makefile for makex on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:15:14 UTC 2021

all: makex makej makek rxscan xtorx



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

MAKEX_LIB = makex_lib.a
MAKEJ_LIB = makej_lib.a
MAKEK_LIB = makek_lib.a
RXSCAN_LIB = rxscan_lib.a
XTORX_LIB = xtorx_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

makex : makex.f $(MAKEX_LIB) $(SLIB) ../includes/dimpar.h ../includes/errflg.h ../includes/global.h ../includes/makex.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) makex.f $(MAKEX_LIB) $(SLIB) -o makex
	rm -f makex.o

makej : makej.f $(MAKEJ_LIB) $(SLIB) ../includes/dimpar.h ../includes/makex.h
	$(FTN) $(FFLAGS) makej.f $(MAKEJ_LIB) $(SLIB) -o makej
	rm -f makej.o

makek : makek.f $(MAKEK_LIB) $(SLIB) ../includes/dimpar.h ../includes/global.h ../includes/makex.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) makek.f $(MAKEK_LIB) $(SLIB) -o makek
	rm -f makek.o

rxscan : rxscan.f $(RXSCAN_LIB) $(SLIB) ../includes/dimpar.h ../includes/makex.h
	$(FTN) $(FFLAGS) rxscan.f $(RXSCAN_LIB) $(SLIB) -o rxscan
	rm -f rxscan.o

xtorx : xtorx.f $(XTORX_LIB) $(SLIB) ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
	$(FTN) $(FFLAGS) xtorx.f $(XTORX_LIB) $(SLIB) -o xtorx
	rm -f xtorx.o

$(MAKEX_LIB) : \
	$(MAKEX_LIB)(closem.o) \
	$(MAKEX_LIB)(getnav.o) \
	$(MAKEX_LIB)(get_rxfiles.o) \
	$(MAKEX_LIB)(getsp3.o) \
	$(MAKEX_LIB)(inqf.o) \
	$(MAKEX_LIB)(inscen.o) \
	$(MAKEX_LIB)(isgood.o) \
	$(MAKEX_LIB)(mversn.o) \
	$(MAKEX_LIB)(openf.o) \
	$(MAKEX_LIB)(rbatch.o) \
	$(MAKEX_LIB)(rcoord.o) \
	$(MAKEX_LIB)(rhead.o) \
	$(MAKEX_LIB)(set_dcb_flag.o) \
	$(MAKEX_LIB)(settim.o) \
	$(MAKEX_LIB)(settyp.o) \
	$(MAKEX_LIB)(stnclk.o) \
	$(MAKEX_LIB)(tagchk.o) \
	$(MAKEX_LIB)(unpad.o) \
	$(MAKEX_LIB)(wclock.o) \
	$(MAKEX_LIB)(wxhead.o)
	ranlib makex_lib.a

$(MAKEJ_LIB) : \
	$(MAKEJ_LIB)(jreadc.o) \
	$(MAKEJ_LIB)(avgclk.o) \
	$(MAKEJ_LIB)(allanv.o) \
	$(MAKEJ_LIB)(clkfit.o) \
	$(MAKEJ_LIB)(get_navclk.o) \
	$(MAKEJ_LIB)(rate_est.o) \
	$(MAKEJ_LIB)(j_from_c.o) \
	$(MAKEJ_LIB)(j_from_nav.o) \
	$(MAKEJ_LIB)(j_from_sp3.o) \
	$(MAKEJ_LIB)(mversn.o) \
	$(MAKEJ_LIB)(pickfns.o)
	ranlib makej_lib.a

$(MAKEK_LIB) : \
	$(MAKEK_LIB)(closem.o) \
	$(MAKEK_LIB)(getnav.o) \
	$(MAKEK_LIB)(getpr.o) \
	$(MAKEK_LIB)(rcoord.o) \
	$(MAKEK_LIB)(stnclk.o) \
	$(MAKEK_LIB)(stncrd.o) \
	$(MAKEK_LIB)(wclock.o) \
	$(MAKEK_LIB)(mversn.o)
	ranlib makek_lib.a

$(RXSCAN_LIB) : \
	$(RXSCAN_LIB)(mversn.o)
	ranlib rxscan_lib.a

$(XTORX_LIB) : \
	$(XTORX_LIB)(mversn.o) \
	$(XTORX_LIB)(scanx.o) \
	$(XTORX_LIB)(settyp.o)
	ranlib xtorx_lib.a

$(MAKEX_LIB)(closem.o): ../includes/makex.h
$(MAKEX_LIB)(getnav.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(get_rxfiles.o): ../includes/makex.h
$(MAKEX_LIB)(getsp3.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(inqf.o): ../includes/makex.h
$(MAKEX_LIB)(inscen.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(isgood.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(openf.o): ../includes/makex.h
$(MAKEX_LIB)(rbatch.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(rcoord.o): ../includes/makex.h
$(MAKEX_LIB)(rhead.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(settim.o): ../includes/makex.h
$(MAKEX_LIB)(settyp.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEX_LIB)(stnclk.o): ../includes/makex.h
$(MAKEX_LIB)(tagchk.o): ../includes/makex.h
$(MAKEX_LIB)(wxhead.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEJ_LIB)(jreadc.o): ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(MAKEJ_LIB)(avgclk.o): ../includes/dimpar.h
$(MAKEJ_LIB)(clkfit.o): ../includes/makex.h
$(MAKEJ_LIB)(get_navclk.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEJ_LIB)(rate_est.o): ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(MAKEJ_LIB)(j_from_c.o): ../includes/dimpar.h ../includes/errflg.h ../includes/makex.h
$(MAKEJ_LIB)(j_from_nav.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEJ_LIB)(j_from_sp3.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEK_LIB)(closem.o): ../includes/makex.h
$(MAKEK_LIB)(getnav.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEK_LIB)(getpr.o): ../includes/dimpar.h ../includes/makex.h
$(MAKEK_LIB)(rcoord.o): ../includes/makex.h
$(MAKEK_LIB)(stnclk.o): ../includes/makex.h
$(XTORX_LIB)(scanx.o): ../includes/dimpar.h ../includes/makex.h
$(XTORX_LIB)(settyp.o): ../includes/dimpar.h ../includes/makex.h

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
