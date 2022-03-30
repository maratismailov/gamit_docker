# Makefile for fica on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:18:37 UTC 2021

all: ficachop ficacomp ficalabl ficamerg ficamove ficapull ficascan ficaslim ficaspan ficaweek fica2rx



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

FICASCAN_LIB = ficascan_lib.a
FICASPAN_LIB = ficaspan_lib.a
FIC2RX_LIB = fic2rx_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

ficachop : ficachop.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficachop.f  $(SLIB) -o ficachop
	rm -f ficachop.o

ficacomp : ficacomp.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficacomp.f  $(SLIB) -o ficacomp
	rm -f ficacomp.o

ficalabl : ficalabl.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficalabl.f  $(SLIB) -o ficalabl
	rm -f ficalabl.o

ficamerg : ficamerg.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficamerg.f  $(SLIB) -o ficamerg
	rm -f ficamerg.o

ficamove : ficamove.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficamove.f  $(SLIB) -o ficamove
	rm -f ficamove.o

ficapull : ficapull.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficapull.f  $(SLIB) -o ficapull
	rm -f ficapull.o

ficascan : ficascan.f $(FICASCAN_LIB) $(SLIB) ../includes/makex.h
	$(FTN) $(FFLAGS) ficascan.f $(FICASCAN_LIB) $(SLIB) -o ficascan
	rm -f ficascan.o

ficaslim : ficaslim.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficaslim.f  $(SLIB) -o ficaslim
	rm -f ficaslim.o

ficaspan : ficaspan.f $(FICASPAN_LIB) $(SLIB) ../includes/makex.h
	$(FTN) $(FFLAGS) ficaspan.f $(FICASPAN_LIB) $(SLIB) -o ficaspan
	rm -f ficaspan.o

ficaweek : ficaweek.f  $(SLIB) 
	$(FTN) $(FFLAGS) ficaweek.f  $(SLIB) -o ficaweek
	rm -f ficaweek.o

fica2rx : fic2rx.f $(FIC2RX_LIB) $(SLIB) ../includes/dimpar.h ../includes/makex.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) fic2rx.f $(FIC2RX_LIB) $(SLIB) -o fic2rx
	rm -f fic2rx.o

$(FICASCAN_LIB) : \
	$(FICASCAN_LIB)(blk101.o)
	ranlib ficascan_lib.a

$(FICASPAN_LIB) : \
	$(FICASPAN_LIB)(dofica.o) \
	$(FICASPAN_LIB)(blk101.o) \
	$(FICASPAN_LIB)(blk1380.o) \
	$(FICASPAN_LIB)(blk6.o) \
	$(FICASPAN_LIB)(blk70.o) \
	$(FICASPAN_LIB)(blk80.o) \
	$(FICASPAN_LIB)(blk9.o) \
	$(FICASPAN_LIB)(blk401.o) \
	$(FICASPAN_LIB)(blk1001.o) \
	$(FICASPAN_LIB)(blk670.o) \
	$(FICASPAN_LIB)(blk1080.o) \
	$(FICASPAN_LIB)(blk1180.o) \
	$(FICASPAN_LIB)(blk1280.o)
	ranlib ficaspan_lib.a

$(FIC2RX_LIB) : \
	$(FIC2RX_LIB)(blk1001.o) \
	$(FIC2RX_LIB)(blk101.o) \
	$(FIC2RX_LIB)(blk1101.o) \
	$(FIC2RX_LIB)(blk1201.o) \
	$(FIC2RX_LIB)(blk1301.o) \
	$(FIC2RX_LIB)(blk1080.o) \
	$(FIC2RX_LIB)(blk1180.o) \
	$(FIC2RX_LIB)(blk1280.o) \
	$(FIC2RX_LIB)(blk1380.o) \
	$(FIC2RX_LIB)(blk401.o) \
	$(FIC2RX_LIB)(blk6.o) \
	$(FIC2RX_LIB)(blk670.o) \
	$(FIC2RX_LIB)(blk70.o) \
	$(FIC2RX_LIB)(blk80.o) \
	$(FIC2RX_LIB)(blk9.o) \
	$(FIC2RX_LIB)(dofica.o) \
	$(FIC2RX_LIB)(rficahd.o)
	ranlib fic2rx_lib.a

$(FICASCAN_LIB)(blk101.o): ../includes/makex.h
$(FICASPAN_LIB)(dofica.o): ../includes/makex.h
$(FICASPAN_LIB)(blk101.o): ../includes/makex.h
$(FICASPAN_LIB)(blk1380.o): ../includes/makex.h
$(FICASPAN_LIB)(blk6.o): ../includes/makex.h
$(FICASPAN_LIB)(blk70.o): ../includes/makex.h
$(FICASPAN_LIB)(blk80.o): ../includes/makex.h
$(FICASPAN_LIB)(blk9.o): ../includes/makex.h
$(FICASPAN_LIB)(blk401.o): ../includes/makex.h
$(FICASPAN_LIB)(blk1001.o): ../includes/makex.h
$(FICASPAN_LIB)(blk670.o): ../includes/makex.h
$(FICASPAN_LIB)(blk1080.o): ../includes/makex.h
$(FICASPAN_LIB)(blk1180.o): ../includes/makex.h
$(FICASPAN_LIB)(blk1280.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1001.o): ../includes/makex.h
$(FIC2RX_LIB)(blk101.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1101.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1201.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1301.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1080.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1180.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1280.o): ../includes/makex.h
$(FIC2RX_LIB)(blk1380.o): ../includes/makex.h
$(FIC2RX_LIB)(blk401.o): ../includes/makex.h
$(FIC2RX_LIB)(blk6.o): ../includes/makex.h
$(FIC2RX_LIB)(blk670.o): ../includes/makex.h
$(FIC2RX_LIB)(blk70.o): ../includes/makex.h
$(FIC2RX_LIB)(blk80.o): ../includes/makex.h
$(FIC2RX_LIB)(blk9.o): ../includes/makex.h
$(FIC2RX_LIB)(dofica.o): ../includes/makex.h
$(FIC2RX_LIB)(rficahd.o): ../includes/makex.h

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
