# Makefile for help on SunOS 5.7 by unimake at Tuesday February 20 17:13:58 EST 2001

all: help_fonda



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

SLIB = ../com/foncom.a 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

help_fonda : help_fonda.f  $(SLIB) 
	$(FTN) $(FFLAGS) help_fonda.f  $(SLIB) -o fonda
	rm -f help_fonda.o


.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
