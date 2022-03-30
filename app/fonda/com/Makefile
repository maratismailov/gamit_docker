# Makefile for com on SunOS 5.7 by unimake at Wednesday May 30 20:31:49 EST 2001

all: foncom



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

FONCOM = foncom.a
SLIB = 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

foncom : $(FONCOM)

$(FONCOM) : \
	$(FONCOM)(ata.o) \
	$(FONCOM)(atwa.o) \
	$(FONCOM)(axb.o) \
	$(FONCOM)(blank.o) \
	$(FONCOM)(cdate.o) \
	$(FONCOM)(chi2_dis.o) \
	$(FONCOM)(clnbuf.o) \
	$(FONCOM)(copy1d.o) \
	$(FONCOM)(count_arg.o) \
	$(FONCOM)(fcheck.o) \
	$(FONCOM)(ferror.o) \
	$(FONCOM)(geotab.o) \
	$(FONCOM)(geoxyz.o) \
	$(FONCOM)(get_path.o) \
	$(FONCOM)(getcmd.o) \
	$(FONCOM)(getdat.o) \
	$(FONCOM)(gettim.o) \
	$(FONCOM)(getusr.o) \
	$(FONCOM)(iclarg.o) \
	$(FONCOM)(isort.o) \
	$(FONCOM)(isort2.o) \
	$(FONCOM)(julday.o) \
	$(FONCOM)(latwa.o) \
	$(FONCOM)(leapyr.o) \
	$(FONCOM)(latwa_indx.o) \
	$(FONCOM)(laxb_indx.o) \
	$(FONCOM)(lift_arg.o) \
	$(FONCOM)(lowerc.o) \
	$(FONCOM)(match_name.o) \
	$(FONCOM)(museum.o) \
	$(FONCOM)(nblen.o) \
	$(FONCOM)(pline.o) \
	$(FONCOM)(remedy_space.o) \
	$(FONCOM)(rtodms.o) \
	$(FONCOM)(satwa.o) \
	$(FONCOM)(sort_int1.o) \
	$(FONCOM)(sort_intr.o) \
	$(FONCOM)(sort_real1.o) \
	$(FONCOM)(sph_ca.o) \
	$(FONCOM)(sphxyz.o) \
	$(FONCOM)(strjst.o) \
	$(FONCOM)(trans.o) \
	$(FONCOM)(upperc.o) \
	$(FONCOM)(wdms.o) \
	$(FONCOM)(zero1d.o) \
	$(FONCOM)(mallocg_sol.o)
	ranlib foncom.a


.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
