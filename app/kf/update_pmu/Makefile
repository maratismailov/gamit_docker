# Makefile for update_pmu on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:40:11 UTC 2021

all: update_pmu update_pmu_tab fix_ut1_tide fix_apr gsfc_to_iris mit_to_iris gamit_to_iris compare_pmu excite_pmu nmc_to_iris usno_to_iris mjd_to_date utc_to_tai igs_to_pmu extend_iers



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

SLIB = ../utils/gen_ssu_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

update_pmu : update_pmu.f  $(SLIB) ../includes/kalman_param.h ../includes/obs_apr.h ../includes/obs_data.h ../includes/obs_header.h
	$(FTN) $(FFLAGS) update_pmu.f  $(SLIB) -o update_pmu
	rm -f update_pmu.o

update_pmu_tab : update_pmu_tab.f  $(SLIB) ../includes/kalman_param.h ../includes/obs_apr.h ../includes/obs_data.h ../includes/obs_header.h
	$(FTN) $(FFLAGS) update_pmu_tab.f  $(SLIB) -o update_pmu_tab
	rm -f update_pmu_tab.o

fix_ut1_tide : fix_ut1_tide.f  $(SLIB) ../includes/kalman_param.h ../includes/obs_data.h ../includes/obs_header.h
	$(FTN) $(FFLAGS) fix_ut1_tide.f  $(SLIB) -o fix_ut1_tide
	rm -f fix_ut1_tide.o

fix_apr : fix_apr.f  $(SLIB) ../includes/kalman_param.h ../includes/obs_header.h
	$(FTN) $(FFLAGS) fix_apr.f  $(SLIB) -o fix_apr
	rm -f fix_apr.o

gsfc_to_iris : gsfc_to_iris.f  $(SLIB) 
	$(FTN) $(FFLAGS) gsfc_to_iris.f  $(SLIB) -o gsfc_to_iris
	rm -f gsfc_to_iris.o

mit_to_iris : mit_to_iris.f  $(SLIB) 
	$(FTN) $(FFLAGS) mit_to_iris.f  $(SLIB) -o mit_to_iris
	rm -f mit_to_iris.o

gamit_to_iris : gamit_to_iris.f  $(SLIB) 
	$(FTN) $(FFLAGS) gamit_to_iris.f  $(SLIB) -o gamit_to_iris
	rm -f gamit_to_iris.o

compare_pmu : compare_pmu.f  $(SLIB) compare_pmu.h
	$(FTN) $(FFLAGS) compare_pmu.f  $(SLIB) -o compare_pmu
	rm -f compare_pmu.o

excite_pmu : excite_pmu.f  $(SLIB) excite_pmu.h
	$(FTN) $(FFLAGS) excite_pmu.f  $(SLIB) -o excite_pmu
	rm -f excite_pmu.o

nmc_to_iris : nmc_to_iris.f  $(SLIB) ../includes/const_param.h
	$(FTN) $(FFLAGS) nmc_to_iris.f  $(SLIB) -o nmc_to_iris
	rm -f nmc_to_iris.o

usno_to_iris : usno_to_iris.f  $(SLIB) 
	$(FTN) $(FFLAGS) usno_to_iris.f  $(SLIB) -o usno_to_iris
	rm -f usno_to_iris.o

mjd_to_date : mjd_to_date.f  $(SLIB) 
	$(FTN) $(FFLAGS) mjd_to_date.f  $(SLIB) -o mjd_to_date
	rm -f mjd_to_date.o

utc_to_tai : utc_to_tai.f  $(SLIB) 
	$(FTN) $(FFLAGS) utc_to_tai.f  $(SLIB) -o utc_to_tai
	rm -f utc_to_tai.o

igs_to_pmu : igs_to_pmu.f  $(SLIB) 
	$(FTN) $(FFLAGS) igs_to_pmu.f  $(SLIB) -o igs_to_pmu
	rm -f igs_to_pmu.o

extend_iers : extend_iers.f  $(SLIB) 
	$(FTN) $(FFLAGS) extend_iers.f  $(SLIB) -o extend_iers
	rm -f extend_iers.o


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
