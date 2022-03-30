# Makefile for maked on SunOS 5.7 by unimake at Wednesday May 30 20:29:20 EST 2001

all: maked



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

MAKED_LIB = maked_lib.a
SLIB = ./maked_lib.a ../com/foncom.a ../kfcomlib/kf_com_lib.a 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

maked : maked.f $(MAKED_LIB) $(SLIB) maked.fti
	$(FTN) $(FFLAGS) maked.f $(MAKED_LIB) $(SLIB) -o maked
	rm -f maked.o

$(MAKED_LIB) : \
	$(MAKED_LIB)(basazm.o) \
	$(MAKED_LIB)(chk_quake_list.o) \
	$(MAKED_LIB)(cmodel.o) \
	$(MAKED_LIB)(compress_dat.o) \
	$(MAKED_LIB)(fonda_out_compress.o) \
	$(MAKED_LIB)(frmobs.o) \
	$(MAKED_LIB)(fmtsft_dat.o) \
	$(MAKED_LIB)(fmtsft_net.o) \
	$(MAKED_LIB)(getcnt.o) \
	$(MAKED_LIB)(getnet.o) \
	$(MAKED_LIB)(grep_data.o) \
	$(MAKED_LIB)(merge_list.o) \
	$(MAKED_LIB)(museum.o) \
	$(MAKED_LIB)(net_head.o) \
	$(MAKED_LIB)(outdat.o) \
	$(MAKED_LIB)(read_bbook_dat.o) \
	$(MAKED_LIB)(read_fondah_dat.o) \
	$(MAKED_LIB)(read_free_net.o) \
	$(MAKED_LIB)(read_gipsy_dat.o) \
	$(MAKED_LIB)(read_glbk_dat.o) \
	$(MAKED_LIB)(read_glbk_full_dat.o) \
	$(MAKED_LIB)(read_glbk_full_net.o) \
	$(MAKED_LIB)(read_hfile_dat.o) \
	$(MAKED_LIB)(read_ipgp_dat.o) \
	$(MAKED_LIB)(read_newbb_dat.o) \
	$(MAKED_LIB)(read_ofile_dat.o) \
	$(MAKED_LIB)(read_ucsd_dat.o) \
	$(MAKED_LIB)(read_usgs_dat.o) \
	$(MAKED_LIB)(read_bbook_net.o) \
	$(MAKED_LIB)(read_glbk_net.o) \
	$(MAKED_LIB)(read_gipsy_net.o) \
	$(MAKED_LIB)(read_ipgp_net.o) \
	$(MAKED_LIB)(read_lfile.o) \
	$(MAKED_LIB)(read_newbb_net.o) \
	$(MAKED_LIB)(read_ref_net.o) \
	$(MAKED_LIB)(read_sinex_dat.o) \
	$(MAKED_LIB)(read_ucsd_net.o) \
	$(MAKED_LIB)(read_usgs_net.o) \
	$(MAKED_LIB)(redinf.o) \
	$(MAKED_LIB)(read_driv.o) \
	$(MAKED_LIB)(unique_nam1.o) \
	$(MAKED_LIB)(unique_name.o) \
	$(MAKED_LIB)(upperc.o) \
	$(MAKED_LIB)(vmodel.o)
	ranlib maked_lib.a

$(MAKED_LIB)(basazm.o): maked.fti
$(MAKED_LIB)(chk_quake_list.o): maked.fti
$(MAKED_LIB)(cmodel.o): maked.fti
$(MAKED_LIB)(compress_dat.o): maked.fti
$(MAKED_LIB)(fonda_out_compress.o): maked.fti
$(MAKED_LIB)(frmobs.o): maked.fti
$(MAKED_LIB)(fmtsft_dat.o): maked.fti
$(MAKED_LIB)(fmtsft_net.o): maked.fti
$(MAKED_LIB)(getcnt.o): maked.fti
$(MAKED_LIB)(getnet.o): maked.fti
$(MAKED_LIB)(net_head.o): maked.fti
$(MAKED_LIB)(outdat.o): maked.fti
$(MAKED_LIB)(read_bbook_dat.o): maked.fti
$(MAKED_LIB)(read_fondah_dat.o): maked.fti
$(MAKED_LIB)(read_free_net.o): maked.fti
$(MAKED_LIB)(read_gipsy_dat.o): maked.fti
$(MAKED_LIB)(read_glbk_dat.o): maked.fti
$(MAKED_LIB)(read_glbk_full_dat.o): maked.fti
$(MAKED_LIB)(read_glbk_full_net.o): maked.fti
$(MAKED_LIB)(read_hfile_dat.o): maked.fti
$(MAKED_LIB)(read_ipgp_dat.o): maked.fti
$(MAKED_LIB)(read_newbb_dat.o): maked.fti
$(MAKED_LIB)(read_ofile_dat.o): maked.fti
$(MAKED_LIB)(read_ucsd_dat.o): maked.fti
$(MAKED_LIB)(read_usgs_dat.o): maked.fti
$(MAKED_LIB)(read_ipgp_net.o): maked.fti
$(MAKED_LIB)(read_ref_net.o): maked.fti
$(MAKED_LIB)(read_sinex_dat.o): ../includes/glb_hdr_def.h ../includes/htoglb_comm.h ../includes/kalman_param.h ../includes/sln_def.h maked.fti
$(MAKED_LIB)(redinf.o): maked.fti
$(MAKED_LIB)(read_driv.o): maked.fti
$(MAKED_LIB)(vmodel.o): maked.fti

.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
