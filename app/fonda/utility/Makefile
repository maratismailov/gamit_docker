# Makefile for utility on SunOS 5.7 by unimake at Tuesday February 20 17:13:42 EST 2001

all: align_frame calc_bsl_azm compare_coor get_line get_scater get_v_rel net_update plot_line rem_model syn_bsl_rat



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

SLIB = ../com/foncom.a 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

align_frame : align_frame.f  $(SLIB) 
	$(FTN) $(FFLAGS) align_frame.f  $(SLIB) -o align_frame
	rm -f align_frame.o

calc_bsl_azm : calc_bsl_azm.f  $(SLIB) 
	$(FTN) $(FFLAGS) calc_bsl_azm.f  $(SLIB) -o calc_bsl_azm
	rm -f calc_bsl_azm.o

compare_coor : compare_coor.f  $(SLIB) 
	$(FTN) $(FFLAGS) compare_coor.f  $(SLIB) -o compare_coor
	rm -f compare_coor.o

get_line : get_line.f  $(SLIB) 
	$(FTN) $(FFLAGS) get_line.f  $(SLIB) -o get_line
	rm -f get_line.o

get_scater : get_scater.f  $(SLIB) 
	$(FTN) $(FFLAGS) get_scater.f  $(SLIB) -o get_scater
	rm -f get_scater.o

get_v_rel : get_v_rel.f  $(SLIB) 
	$(FTN) $(FFLAGS) get_v_rel.f  $(SLIB) -o get_v_rel
	rm -f get_v_rel.o

net_update : net_update.f  $(SLIB) 
	$(FTN) $(FFLAGS) net_update.f  $(SLIB) -o net_update
	rm -f net_update.o

plot_line : plot_line.f  $(SLIB) 
	$(FTN) $(FFLAGS) plot_line.f  $(SLIB) -o plot_line
	rm -f plot_line.o

rem_model : rem_model.f  $(SLIB) 
	$(FTN) $(FFLAGS) rem_model.f  $(SLIB) -o rem_model
	rm -f rem_model.o

syn_bsl_rat : syn_bsl_rat.f  $(SLIB) 
	$(FTN) $(FFLAGS) syn_bsl_rat.f  $(SLIB) -o syn_bsl_rat
	rm -f syn_bsl_rat.o


.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
