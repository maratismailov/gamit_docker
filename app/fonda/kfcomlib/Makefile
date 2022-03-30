# Makefile for kfcomlib on SunOS 5.7 by unimake at Thursday February 22 13:59:29 EST 2001

all: kf_com_lib



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

KF_COM_LIB = kf_com_lib.a
SLIB = 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

kf_com_lib : $(KF_COM_LIB)

$(KF_COM_LIB) : \
	$(KF_COM_LIB)(casefold.o) \
	$(KF_COM_LIB)(trimlead.o) \
	$(KF_COM_LIB)(yds_to_jd.o) \
	$(KF_COM_LIB)(trim_lead.o) \
	$(KF_COM_LIB)(report_error.o) \
	$(KF_COM_LIB)(invert_vis.o) \
	$(KF_COM_LIB)(jacobi.o) \
	$(KF_COM_LIB)(get_command.o) \
	$(KF_COM_LIB)(trimlead.o) \
	$(KF_COM_LIB)(dwmov.o) \
	$(KF_COM_LIB)(dwvmv.o) \
	$(KF_COM_LIB)(short_lod.o) \
	$(KF_COM_LIB)(gauss_elim.o) \
	$(KF_COM_LIB)(dwsmy.o) \
	$(KF_COM_LIB)(dwsmy.o) \
	$(KF_COM_LIB)(dwmab.o) \
	$(KF_COM_LIB)(dwswp.o) \
	$(KF_COM_LIB)(dwpiv.o) \
	$(KF_COM_LIB)(pmu_main_par.o) \
	$(KF_COM_LIB)(scale_matrix.o) \
	$(KF_COM_LIB)(trimlen.o) \
	$(KF_COM_LIB)(ymdhms_to_jd.o) \
	$(KF_COM_LIB)(decode_snxb.o) \
	$(KF_COM_LIB)(decsnx_site.o) \
	$(KF_COM_LIB)(decsnx_soln.o) \
	$(KF_COM_LIB)(preread_snx.o) \
	$(KF_COM_LIB)(snx_finbl.o) \
	$(KF_COM_LIB)(decsnx_file.o) \
	$(KF_COM_LIB)(decsnx_input.o) \
	$(KF_COM_LIB)(GetWord.o) \
	$(KF_COM_LIB)(report_stat.o) \
	$(KF_COM_LIB)(mdian2.o) \
	$(KF_COM_LIB)(check_ascii.o) \
	$(KF_COM_LIB)(caseunfold.o) \
	$(KF_COM_LIB)(systime.o) \
	$(KF_COM_LIB)(jd_to_decyrs.o) \
	$(KF_COM_LIB)(jd_to_ymdhms.o)
	ranlib kf_com_lib.a

$(KF_COM_LIB)(pmu_main_par.o): ../includes/const_param.h
$(KF_COM_LIB)(decode_snxb.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h
$(KF_COM_LIB)(decsnx_site.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h
$(KF_COM_LIB)(decsnx_soln.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h
$(KF_COM_LIB)(preread_snx.o): ../includes/htoglb_comm.h ../includes/kalman_param.h
$(KF_COM_LIB)(decsnx_file.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h
$(KF_COM_LIB)(decsnx_input.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h

.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
