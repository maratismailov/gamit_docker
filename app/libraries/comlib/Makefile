# Makefile for comlib on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:12:04 UTC 2021

all: com_lib



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

COM_LIB = com_lib.a
SLIB = 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

com_lib : $(COM_LIB)

$(COM_LIB) : \
	$(COM_LIB)(sd_compbd.o) \
	$(COM_LIB)(GetWord.o) \
	$(COM_LIB)(MoveChars.o) \
	$(COM_LIB)(addressof.o) \
	$(COM_LIB)(bad_option.o) \
	$(COM_LIB)(casefold.o) \
	$(COM_LIB)(caseunfold.o) \
	$(COM_LIB)(check_ascii.o) \
	$(COM_LIB)(check_reset.o) \
	$(COM_LIB)(clearbuffer.o) \
	$(COM_LIB)(cversn.o) \
	$(COM_LIB)(decimaltoint.o) \
	$(COM_LIB)(decyrs_to_jd.o) \
	$(COM_LIB)(decyrs_to_ydhms.o) \
	$(COM_LIB)(ds2hms.o) \
	$(COM_LIB)(eposn.o) \
	$(COM_LIB)(fmpclose.o) \
	$(COM_LIB)(fmpendmask.o) \
	$(COM_LIB)(fmpinitmask.o) \
	$(COM_LIB)(fmpnextmask.o) \
	$(COM_LIB)(fmppurge.o) \
	$(COM_LIB)(fullfilename.o) \
	$(COM_LIB)(gipson.o) \
	$(COM_LIB)(gen_help_nm.o) \
	$(COM_LIB)(get_cmd.o) \
	$(COM_LIB)(get_value.o) \
	$(COM_LIB)(iabs2.o) \
	$(COM_LIB)(ifbrk.o) \
	$(COM_LIB)(jd_to_decyrs.o) \
	$(COM_LIB)(jd_to_yds.o) \
	$(COM_LIB)(jd_to_ymdhms.o) \
	$(COM_LIB)(leapyr.o) \
	$(COM_LIB)(leq.o) \
	$(COM_LIB)(ljust.o) \
	$(COM_LIB)(loglu.o) \
	$(COM_LIB)(lunit.o) \
	$(COM_LIB)(min02.o) \
	$(COM_LIB)(max02.o) \
	$(COM_LIB)(mjd_to_ymdhm.o) \
	$(COM_LIB)(movewords.o) \
	$(COM_LIB)(multiread.o) \
	$(COM_LIB)(null_termn.o) \
	$(COM_LIB)(open_lu.o) \
	$(COM_LIB)(pmut1_oceans.o) \
	$(COM_LIB)(proper_runst.o) \
	$(COM_LIB)(null_termn.o) \
	$(COM_LIB)(open_lu.o) \
	$(COM_LIB)(proper_runst.o) \
	$(COM_LIB)(ray.o) \
	$(COM_LIB)(rcpar.o) \
	$(COM_LIB)(rcpar4.o) \
	$(COM_LIB)(read_line.o) \
	$(COM_LIB)(readd.o) \
	$(COM_LIB)(readd8.o) \
	$(COM_LIB)(report_error.o) \
	$(COM_LIB)(report_stat.o) \
	$(COM_LIB)(sd_comp.o) \
	$(COM_LIB)(sub_char.o) \
	$(COM_LIB)(swap_bytes.o) \
	$(COM_LIB)(tide_angles.o) \
	$(COM_LIB)(timr.o) \
	$(COM_LIB)(trimlead.o) \
	$(COM_LIB)(trimlen.o) \
	$(COM_LIB)(underscore.o) \
	$(COM_LIB)(wline.o) \
	$(COM_LIB)(wmove.o) \
	$(COM_LIB)(writd.o) \
	$(COM_LIB)(writd8.o) \
	$(COM_LIB)(yds_to_jd.o) \
	$(COM_LIB)(ymdhms_to_jd.o) \
	$(COM_LIB)(ymdhms_to_mj.o) \
	$(COM_LIB)(ymd_to_doy.o) \
	$(COM_LIB)(jd_to_decyrs.o) \
	$(COM_LIB)(jd_to_gpst.o) \
	$(COM_LIB)(memassign.o) \
	$(COM_LIB)(memassign8.o) \
	$(COM_LIB)(GEOD_to_XYZ.o) \
	$(COM_LIB)(XYZ_to_GEOD.o) \
	$(COM_LIB)(XYZ_to_NEU.o) \
	$(COM_LIB)(mean_pole.o) \
	$(COM_LIB)(name_to_blk.o) \
	$(COM_LIB)(bit_util_gftn.o) \
	$(COM_LIB)(cand_gftn.o) \
	$(COM_LIB)(cor_gftn.o) \
	$(COM_LIB)(fmpopen_gftn.o) \
	$(COM_LIB)(getdir_gftn.o) \
	$(COM_LIB)(getmac_gftn.o) \
	$(COM_LIB)(inkey_gftn.o) \
	$(COM_LIB)(irename_gftn.o) \
	$(COM_LIB)(pickfn_gftn.o) \
	$(COM_LIB)(fseekg_gftn.o) \
	$(COM_LIB)(fmprename_gftn.o) \
	$(COM_LIB)(fmpshortname_gftn.o) \
	$(COM_LIB)(fmprunprog_gftn.o) \
	$(COM_LIB)(ferror_gftn.o) \
	$(COM_LIB)(systime_gftn.o) \
	$(COM_LIB)(execute_gftn.o)
	ranlib com_lib.a

$(COM_LIB)(GEOD_to_XYZ.o): ../includes/const_param.h
$(COM_LIB)(XYZ_to_GEOD.o): ../includes/const_param.h

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
