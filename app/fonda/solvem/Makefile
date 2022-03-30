# Makefile for solvem on SunOS 5.7 by unimake at Tuesday February 20 17:13:34 EST 2001

all: solvem



# ASSIGMENTS

FFLAGS = -O3 -u -Nl100
FTN = f77 
CFLAGS = -g -c
CC = cc

SOLVEM_ALL = solvem_all.a
SLIB = ../com/foncom.a 
X11LIB = /usr/openwin/lib
II = /usr/openwin/share/include

solvem : solvem.f $(SOLVEM_ALL) $(SLIB) solvem.fti
	$(FTN) $(FFLAGS) solvem.f $(SOLVEM_ALL) $(SLIB) -o solvem
	rm -f solvem.o

$(SOLVEM_ALL) : \
	$(SOLVEM_ALL)(addaux.o) \
	$(SOLVEM_ALL)(add_quake_c.o) \
	$(SOLVEM_ALL)(azimuc.o) \
	$(SOLVEM_ALL)(basazm.o) \
	$(SOLVEM_ALL)(baslc.o) \
	$(SOLVEM_ALL)(calc_omc_block.o) \
	$(SOLVEM_ALL)(calc_omc_full.o) \
	$(SOLVEM_ALL)(calc_res.o) \
	$(SOLVEM_ALL)(calstr.o) \
	$(SOLVEM_ALL)(centroid.o) \
	$(SOLVEM_ALL)(chk_epi_list.o) \
	$(SOLVEM_ALL)(chk_quake_list.o) \
	$(SOLVEM_ALL)(cholsk.o) \
	$(SOLVEM_ALL)(chknrm.o) \
	$(SOLVEM_ALL)(cpsoln.o) \
	$(SOLVEM_ALL)(dotd.o) \
	$(SOLVEM_ALL)(delaunay.o) \
	$(SOLVEM_ALL)(exert_link.o) \
	$(SOLVEM_ALL)(fill_apr_mtx.o) \
	$(SOLVEM_ALL)(filln.o) \
	$(SOLVEM_ALL)(filln_block.o) \
	$(SOLVEM_ALL)(filln_full.o) \
	$(SOLVEM_ALL)(frmnor.o) \
	$(SOLVEM_ALL)(frmnor_seqe.o) \
	$(SOLVEM_ALL)(geocnt.o) \
	$(SOLVEM_ALL)(geovcv.o) \
	$(SOLVEM_ALL)(get_iom_list.o) \
	$(SOLVEM_ALL)(getjac.o) \
	$(SOLVEM_ALL)(getomc.o) \
	$(SOLVEM_ALL)(getsit.o) \
	$(SOLVEM_ALL)(getpri.o) \
	$(SOLVEM_ALL)(get_qk_adj.o) \
	$(SOLVEM_ALL)(get_quake_list.o) \
	$(SOLVEM_ALL)(get_refractn.o) \
	$(SOLVEM_ALL)(get_res_full.o) \
	$(SOLVEM_ALL)(get_res_sgl.o) \
	$(SOLVEM_ALL)(gmt_sketch.o) \
	$(SOLVEM_ALL)(gmt_file.o) \
	$(SOLVEM_ALL)(helinv.o) \
	$(SOLVEM_ALL)(hwrite.o) \
	$(SOLVEM_ALL)(inner_sln.o) \
	$(SOLVEM_ALL)(initn.o) \
	$(SOLVEM_ALL)(lsq_soln.o) \
	$(SOLVEM_ALL)(lswght.o) \
	$(SOLVEM_ALL)(model_sln.o) \
	$(SOLVEM_ALL)(museum.o) \
	$(SOLVEM_ALL)(norms.o) \
	$(SOLVEM_ALL)(nrmscl.o) \
	$(SOLVEM_ALL)(obsred.o) \
	$(SOLVEM_ALL)(outer_sln.o) \
	$(SOLVEM_ALL)(output.o) \
	$(SOLVEM_ALL)(poscor.o) \
	$(SOLVEM_ALL)(prechk.o) \
	$(SOLVEM_ALL)(prechk_full.o) \
	$(SOLVEM_ALL)(propcv.o) \
	$(SOLVEM_ALL)(put_apr_wght.o) \
	$(SOLVEM_ALL)(put_fix.o) \
	$(SOLVEM_ALL)(put_link.o) \
	$(SOLVEM_ALL)(put_assign.o) \
	$(SOLVEM_ALL)(read_obs_block.o) \
	$(SOLVEM_ALL)(read_obs_full.o) \
	$(SOLVEM_ALL)(readrv.o) \
	$(SOLVEM_ALL)(recod.o) \
	$(SOLVEM_ALL)(reduce_para.o) \
	$(SOLVEM_ALL)(remedy.o) \
	$(SOLVEM_ALL)(reordr_indx.o) \
	$(SOLVEM_ALL)(residl.o) \
	$(SOLVEM_ALL)(res_list.o) \
	$(SOLVEM_ALL)(rotsln.o) \
	$(SOLVEM_ALL)(seqsln.o) \
	$(SOLVEM_ALL)(setlbl.o) \
	$(SOLVEM_ALL)(sln_update.o) \
	$(SOLVEM_ALL)(strain.o) \
	$(SOLVEM_ALL)(submtx.o) \
	$(SOLVEM_ALL)(subtra.o) \
	$(SOLVEM_ALL)(trim_matrix.o) \
	$(SOLVEM_ALL)(update_pointer.o) \
	$(SOLVEM_ALL)(wthead.o)
	ranlib solvem_all.a

$(SOLVEM_ALL)(addaux.o): solvem.fti
$(SOLVEM_ALL)(add_quake_c.o): solvem.fti
$(SOLVEM_ALL)(azimuc.o): solvem.fti
$(SOLVEM_ALL)(basazm.o): solvem.fti
$(SOLVEM_ALL)(baslc.o): solvem.fti
$(SOLVEM_ALL)(calc_omc_block.o): solvem.fti
$(SOLVEM_ALL)(calc_omc_full.o): solvem.fti
$(SOLVEM_ALL)(calc_res.o): solvem.fti
$(SOLVEM_ALL)(calstr.o): solvem.fti
$(SOLVEM_ALL)(chk_epi_list.o): solvem.fti
$(SOLVEM_ALL)(chk_quake_list.o): solvem.fti
$(SOLVEM_ALL)(cpsoln.o): solvem.fti
$(SOLVEM_ALL)(delaunay.o): solvem.fti
$(SOLVEM_ALL)(exert_link.o): solvem.fti
$(SOLVEM_ALL)(filln.o): solvem.fti
$(SOLVEM_ALL)(filln_block.o): solvem.fti
$(SOLVEM_ALL)(filln_full.o): solvem.fti
$(SOLVEM_ALL)(frmnor.o): solvem.fti
$(SOLVEM_ALL)(frmnor_seqe.o): solvem.fti
$(SOLVEM_ALL)(geocnt.o): solvem.fti
$(SOLVEM_ALL)(geovcv.o): solvem.fti
$(SOLVEM_ALL)(get_iom_list.o): solvem.fti
$(SOLVEM_ALL)(getjac.o): solvem.fti
$(SOLVEM_ALL)(getomc.o): solvem.fti
$(SOLVEM_ALL)(getsit.o): solvem.fti
$(SOLVEM_ALL)(getpri.o): solvem.fti
$(SOLVEM_ALL)(get_qk_adj.o): solvem.fti
$(SOLVEM_ALL)(get_quake_list.o): solvem.fti
$(SOLVEM_ALL)(get_refractn.o): solvem.fti
$(SOLVEM_ALL)(get_res_full.o): solvem.fti
$(SOLVEM_ALL)(get_res_sgl.o): solvem.fti
$(SOLVEM_ALL)(gmt_sketch.o): solvem.fti
$(SOLVEM_ALL)(gmt_file.o): solvem.fti
$(SOLVEM_ALL)(hwrite.o): solvem.fti
$(SOLVEM_ALL)(inner_sln.o): solvem.fti
$(SOLVEM_ALL)(initn.o): solvem.fti
$(SOLVEM_ALL)(lsq_soln.o): solvem.fti
$(SOLVEM_ALL)(lswght.o): solvem.fti
$(SOLVEM_ALL)(model_sln.o): solvem.fti
$(SOLVEM_ALL)(norms.o): solvem.fti
$(SOLVEM_ALL)(obsred.o): solvem.fti
$(SOLVEM_ALL)(outer_sln.o): solvem.fti
$(SOLVEM_ALL)(output.o): solvem.fti
$(SOLVEM_ALL)(poscor.o): solvem.fti
$(SOLVEM_ALL)(prechk.o): solvem.fti
$(SOLVEM_ALL)(prechk_full.o): solvem.fti
$(SOLVEM_ALL)(put_apr_wght.o): solvem.fti
$(SOLVEM_ALL)(put_fix.o): solvem.fti
$(SOLVEM_ALL)(put_link.o): solvem.fti
$(SOLVEM_ALL)(put_assign.o): solvem.fti
$(SOLVEM_ALL)(read_obs_block.o): solvem.fti
$(SOLVEM_ALL)(read_obs_full.o): solvem.fti
$(SOLVEM_ALL)(readrv.o): solvem.fti
$(SOLVEM_ALL)(recod.o): solvem.fti
$(SOLVEM_ALL)(reduce_para.o): solvem.fti
$(SOLVEM_ALL)(remedy.o): solvem.fti
$(SOLVEM_ALL)(reordr_indx.o): solvem.fti
$(SOLVEM_ALL)(residl.o): solvem.fti
$(SOLVEM_ALL)(res_list.o): solvem.fti
$(SOLVEM_ALL)(rotsln.o): solvem.fti
$(SOLVEM_ALL)(setlbl.o): solvem.fti
$(SOLVEM_ALL)(strain.o): solvem.fti
$(SOLVEM_ALL)(subtra.o): solvem.fti
$(SOLVEM_ALL)(wthead.o): solvem.fti

.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
