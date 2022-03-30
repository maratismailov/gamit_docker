# Makefile for grdtab on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:19:58 UTC 2021

all: grdtab daf2da dafexpand nao2da atmtoasc atmtide update_atmlg update_vmfg vmftoasc vmfasc2bin atmfilt atl_asc2bin atltoasc toc2da



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

GRDTAB_LIB = grdtab_lib.a
DAF2DA_LIB = daf2da_lib.a
NAO2DA_LIB = nao2da_lib.a
ATMTOASC_LIB = atmtoasc_lib.a
UPDATE_ATMLG_LIB = update_atmlg_lib.a
UPDATE_VMFG_LIB = update_vmfg_lib.a
VMFTOASC_LIB = vmftoasc_lib.a
ATMFILT_LIB = atmfilt_lib.a
ATL_ASC2BIN_LIB = atl_asc2bin_lib.a
ATLTOASC_LIB = atltoasc_lib.a
TOC2DA_LIB = toc2da_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

grdtab : grdtab.f $(GRDTAB_LIB) $(SLIB) ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h ../includes/units.h
	$(FTN) $(FFLAGS) grdtab.f $(GRDTAB_LIB) $(SLIB) -o grdtab
	rm -f grdtab.o

daf2da : daf2da.f $(DAF2DA_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) daf2da.f $(DAF2DA_LIB) $(SLIB) -o daf2da
	rm -f daf2da.o

dafexpand : dafexpand.f  $(SLIB) 
	$(FTN) $(FFLAGS) dafexpand.f  $(SLIB) -o dafexpand
	rm -f dafexpand.o

nao2da : nao2da.f $(NAO2DA_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) nao2da.f $(NAO2DA_LIB) $(SLIB) -o nao2da
	rm -f nao2da.o

atmtoasc : atmtoasc.f $(ATMTOASC_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) atmtoasc.f $(ATMTOASC_LIB) $(SLIB) -o atmtoasc
	rm -f atmtoasc.o

atmtide : atmtide.f  $(SLIB) 
	$(FTN) $(FFLAGS) atmtide.f  $(SLIB) -o atmtide
	rm -f atmtide.o

update_atmlg : update_atmlg.f $(UPDATE_ATMLG_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) update_atmlg.f $(UPDATE_ATMLG_LIB) $(SLIB) -o update_atmlg
	rm -f update_atmlg.o

update_vmfg : update_vmfg.f $(UPDATE_VMFG_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) update_vmfg.f $(UPDATE_VMFG_LIB) $(SLIB) -o update_vmfg
	rm -f update_vmfg.o

vmftoasc : vmftoasc.f $(VMFTOASC_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) vmftoasc.f $(VMFTOASC_LIB) $(SLIB) -o vmftoasc
	rm -f vmftoasc.o

vmfasc2bin : vmfasc2bin.f  $(SLIB) 
	$(FTN) $(FFLAGS) vmfasc2bin.f  $(SLIB) -o vmfasc2bin
	rm -f vmfasc2bin.o

atmfilt : atmfilt.f $(ATMFILT_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) atmfilt.f $(ATMFILT_LIB) $(SLIB) -o atmfilt
	rm -f atmfilt.o

atl_asc2bin : atl_asc2bin.f $(ATL_ASC2BIN_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) atl_asc2bin.f $(ATL_ASC2BIN_LIB) $(SLIB) -o atl_asc2bin
	rm -f atl_asc2bin.o

atltoasc : atltoasc.f $(ATLTOASC_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) atltoasc.f $(ATLTOASC_LIB) $(SLIB) -o atltoasc
	rm -f atltoasc.o

toc2da : toc2da.f $(TOC2DA_LIB) $(SLIB) 
	$(FTN) $(FFLAGS) toc2da.f $(TOC2DA_LIB) $(SLIB) -o toc2da
	rm -f toc2da.o

$(GRDTAB_LIB) : \
	$(GRDTAB_LIB)(grdtab.o) \
	$(GRDTAB_LIB)(grdtab_bd.o) \
	$(GRDTAB_LIB)(gversn.o) \
	$(GRDTAB_LIB)(bilin4.o) \
	$(GRDTAB_LIB)(geoc_to_geod.o) \
	$(GRDTAB_LIB)(geod_to_xyz.o) \
	$(GRDTAB_LIB)(get_atl_list.o) \
	$(GRDTAB_LIB)(get_map_grid.o) \
	$(GRDTAB_LIB)(get_met_list.o) \
	$(GRDTAB_LIB)(get_atl_grid.o) \
	$(GRDTAB_LIB)(get_atml_list.o) \
	$(GRDTAB_LIB)(get_met_grid.o) \
	$(GRDTAB_LIB)(get_otl_list.o) \
	$(GRDTAB_LIB)(get_atml_grid.o) \
	$(GRDTAB_LIB)(get_map_list.o) \
	$(GRDTAB_LIB)(get_otl_grid.o) \
	$(GRDTAB_LIB)(grid_time.o) \
	$(GRDTAB_LIB)(harpos_head.o) \
	$(GRDTAB_LIB)(interp_months.o) \
	$(GRDTAB_LIB)(rd_atl_grid.o) \
	$(GRDTAB_LIB)(rd_atml_list.o) \
	$(GRDTAB_LIB)(rd_atml_grid.o) \
	$(GRDTAB_LIB)(rd_map_list.o) \
	$(GRDTAB_LIB)(rd_otl_grid.o) \
	$(GRDTAB_LIB)(rd_met_grid.o) \
	$(GRDTAB_LIB)(rd_otl_list.o) \
	$(GRDTAB_LIB)(rd_atl_list.o) \
	$(GRDTAB_LIB)(rd_map_grid.o) \
	$(GRDTAB_LIB)(rd_met_list.o) \
	$(GRDTAB_LIB)(reduce_zhd.o) \
	$(GRDTAB_LIB)(xyz_to_geod.o) \
	$(GRDTAB_LIB)(rep_S1S2CoM.o)
	ranlib grdtab_lib.a

$(DAF2DA_LIB) : \
	$(DAF2DA_LIB)(gversn.o)
	ranlib daf2da_lib.a

$(NAO2DA_LIB) : \
	$(NAO2DA_LIB)(gversn.o)
	ranlib nao2da_lib.a

$(ATMTOASC_LIB) : \
	$(ATMTOASC_LIB)(gversn.o) \
	$(ATMTOASC_LIB)(grid_time.o)
	ranlib atmtoasc_lib.a

$(UPDATE_ATMLG_LIB) : \
	$(UPDATE_ATMLG_LIB)(gversn.o) \
	$(UPDATE_ATMLG_LIB)(grid_time.o)
	ranlib update_atmlg_lib.a

$(UPDATE_VMFG_LIB) : \
	$(UPDATE_VMFG_LIB)(gversn.o) \
	$(UPDATE_VMFG_LIB)(grid_time.o)
	ranlib update_vmfg_lib.a

$(VMFTOASC_LIB) : \
	$(VMFTOASC_LIB)(gversn.o) \
	$(VMFTOASC_LIB)(grid_time.o)
	ranlib vmftoasc_lib.a

$(ATMFILT_LIB) : \
	$(ATMFILT_LIB)(gversn.o) \
	$(ATMFILT_LIB)(atm_butter.o) \
	$(ATMFILT_LIB)(filter_atml.o) \
	$(ATMFILT_LIB)(grid_time.o) \
	$(ATMFILT_LIB)(read_atml_headers.o)
	ranlib atmfilt_lib.a

$(ATL_ASC2BIN_LIB) : \
	$(ATL_ASC2BIN_LIB)(gversn.o)
	ranlib atl_asc2bin_lib.a

$(ATLTOASC_LIB) : \
	$(ATLTOASC_LIB)(gversn.o) \
	$(ATLTOASC_LIB)(interp_months.o)
	ranlib atltoasc_lib.a

$(TOC2DA_LIB) : \
	$(TOC2DA_LIB)(gversn.o)
	ranlib toc2da_lib.a

$(GRDTAB_LIB)(grdtab.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h ../includes/units.h
$(GRDTAB_LIB)(grdtab_bd.o): ../includes/grdtab.h
$(GRDTAB_LIB)(get_atl_list.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_map_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_met_list.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_atl_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_atml_list.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_met_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_otl_list.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_atml_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_map_list.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(get_otl_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(rd_atl_grid.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_atml_list.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_atml_grid.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_map_list.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_otl_grid.o): ../includes/dimpar.h ../includes/grdtab.h ../includes/model.h
$(GRDTAB_LIB)(rd_met_grid.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_otl_list.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_map_grid.o): ../includes/grdtab.h
$(GRDTAB_LIB)(rd_met_list.o): ../includes/grdtab.h

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
