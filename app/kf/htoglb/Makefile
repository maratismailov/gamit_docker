# Makefile for htoglb on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:34:33 UTC 2021

all: htoglb_lib htoglb glbtosnx hfupd mstinf chkstinf eqfupd



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

HTOGLB_LIB = htoglb_lib.a
SLIB = ./htoglb_lib.a ../Khandlers/KalObs_hndl_lib.a ../gen_util/gen_util_lib.a ../../gamit/lib//gamit_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

htoglb_lib : $(HTOGLB_LIB)

htoglb : htoglb.f  $(SLIB) ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
	$(FTN) $(FFLAGS) htoglb.f  $(SLIB) -o htoglb
	rm -f htoglb.o

glbtosnx : glbtosnx.f  $(SLIB) ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h htoglb_comm.h
	$(FTN) $(FFLAGS) glbtosnx.f  $(SLIB) -o glbtosnx
	rm -f glbtosnx.o

hfupd : hfupd.f  $(SLIB) ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h hfupd.h htoglb_comm.h
	$(FTN) $(FFLAGS) hfupd.f  $(SLIB) -o hfupd
	rm -f hfupd.o

mstinf : mstinf.f  $(SLIB) ../../gamit/includes/makex.h mstinf.h
	$(FTN) $(FFLAGS) mstinf.f  $(SLIB) -o mstinf
	rm -f mstinf.o

chkstinf : chkstinf.f  $(SLIB) 
	$(FTN) $(FFLAGS) chkstinf.f  $(SLIB) -o chkstinf
	rm -f chkstinf.o

eqfupd : eqfupd.f  $(SLIB) 
	$(FTN) $(FFLAGS) eqfupd.f  $(SLIB) -o eqfupd
	rm -f eqfupd.o

$(HTOGLB_LIB) : \
	$(HTOGLB_LIB)(htog_ut.o) \
	$(HTOGLB_LIB)(glbsx_ut.o) \
	$(HTOGLB_LIB)(hfupd_ut.o) \
	$(HTOGLB_LIB)(hw_util.o) \
	$(HTOGLB_LIB)(hread_util.o) \
	$(HTOGLB_LIB)(init_qparn.o) \
	$(HTOGLB_LIB)(init_qapr_co.o) \
	$(HTOGLB_LIB)(write_ephem.o) \
	$(HTOGLB_LIB)(get_indx_q.o) \
	$(HTOGLB_LIB)(compress_svs.o) \
	$(HTOGLB_LIB)(compress_sit.o) \
	$(HTOGLB_LIB)(read_orient.o) \
	$(HTOGLB_LIB)(read_apr_cov.o) \
	$(HTOGLB_LIB)(skip_htoh.o) \
	$(HTOGLB_LIB)(get_gps_hext.o) \
	$(HTOGLB_LIB)(get_svs_inf.o) \
	$(HTOGLB_LIB)(gsfc_vlbi.o) \
	$(HTOGLB_LIB)(ut_slr.o) \
	$(HTOGLB_LIB)(make_ctog.o) \
	$(HTOGLB_LIB)(decsnx_file.o) \
	$(HTOGLB_LIB)(decsnx_input.o) \
	$(HTOGLB_LIB)(decsnx_site.o) \
	$(HTOGLB_LIB)(decsnx_sat.o) \
	$(HTOGLB_LIB)(decsnx_soln.o) \
	$(HTOGLB_LIB)(decode_snxb.o) \
	$(HTOGLB_LIB)(make_snx.o) \
	$(HTOGLB_LIB)(preread_snx.o) \
	$(HTOGLB_LIB)(snx_finbl.o) \
	$(HTOGLB_LIB)(create_snx.o) \
	$(HTOGLB_LIB)(cresnx_file.o) \
	$(HTOGLB_LIB)(cresnx_input.o) \
	$(HTOGLB_LIB)(cresnx_site.o) \
	$(HTOGLB_LIB)(cresnx_sat.o) \
	$(HTOGLB_LIB)(cresnx_soln.o) \
	$(HTOGLB_LIB)(gen_snxname.o) \
	$(HTOGLB_LIB)(jacobi.o) \
	$(HTOGLB_LIB)(SMod_to_full.o) \
	$(HTOGLB_LIB)(remove_neq.o) \
	$(HTOGLB_LIB)(pack_neq.o)
	ranlib htoglb_lib.a

$(HTOGLB_LIB)(htog_ut.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(glbsx_ut.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h ../includes/svinf_def.h htoglb_comm.h
$(HTOGLB_LIB)(hfupd_ut.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h hfupd.h htoglb_comm.h
$(HTOGLB_LIB)(hw_util.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sinf_def.h ../includes/sln_def.h ../includes/svinf_def.h htoglb_comm.h
$(HTOGLB_LIB)(hread_util.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(init_qparn.o): ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(init_qapr_co.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(write_ephem.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(compress_svs.o): ../includes/const_param.h ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(compress_sit.o): ../includes/const_param.h ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(read_orient.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(read_apr_cov.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(get_svs_inf.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(gsfc_vlbi.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(ut_slr.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(make_ctog.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(decsnx_file.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(decsnx_input.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(decsnx_site.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(decsnx_sat.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h ../includes/svinf_def.h htoglb_comm.h
$(HTOGLB_LIB)(decsnx_soln.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(decode_snxb.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(make_snx.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(preread_snx.o): ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(create_snx.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(cresnx_file.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(cresnx_input.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(cresnx_site.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(cresnx_sat.o): ../includes/const_param.h ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/svinf_def.h htoglb_comm.h
$(HTOGLB_LIB)(cresnx_soln.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(gen_snxname.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(SMod_to_full.o): ../includes/const_param.h ../includes/kalman_param.h htoglb_comm.h
$(HTOGLB_LIB)(remove_neq.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h
$(HTOGLB_LIB)(pack_neq.o): ../includes/glb_hdr_def.h ../includes/kalman_param.h ../includes/sln_def.h htoglb_comm.h

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
