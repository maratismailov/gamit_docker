# Makefile for trackRT on Linux 2.6.31-20-generic by unimake at Mon Apr 19 11:11:15 EDT 2010

all: trackRT_lib trackRTr trackRT trackRTB test_trackRT

### GCC V4 gfortran flags ###
FTN = gfortran
FFLAGS = -O3 -Wuninitialized -fno-f2c -ffast-math -fno-automatic -fno-backslash
CC = gcc
CFLAGS = -g  
#### For 64-bit machines:
#CFLAGS = -g -DADD64BIT

TRACKRT_LIB = trackRT_lib.a
SLIB = ./trackRT_lib.a ../gen_util/gen_util_lib.a ../../libraries/matrix/kinv_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib
II = /usr/include/X11
QTLIB = /usr/share/qt4/lib
QTINC = /usr/include/qt4

trackRT_lib : $(TRACKRT_LIB)

trackRTr : trackRTr.f  $(SLIB) trackRT.h trackRTObs.h
	$(FTN) $(FFLAGS) trackRTr.f  $(SLIB) -o trackRTr
	rm -f trackRTr.o

trackRT : $(SLIB) trackRT.h trackRTObs.h
	g++ -c -pipe -O2 -Wall -W -D_REENTRANT -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_NETWORK_LIB -DQT_CORE_LIB -DQT_SHARED \
                     -I${QTINC}/mkspecs/linux-g++ -I${QTINC}/QtCore -I${QTINC}/QtNetwork  -I${QTINC} -c trackRTComm.cpp
	$(FTN) $(FFLAGS) trackRTComm.o -L/usr/local/lib -L${QTLIB} -lQtNetwork -lQtCore -lpthread  -lstdc++ $(SLIB) -o trackRT
	rm -f trackRTComm.o

test_trackRT : $(SLIB) trackRT.h trackRTObs.h
	g++ -c -pipe -O2 -Wall -W -D_REENTRANT -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_NETWORK_LIB -DQT_CORE_LIB -DQT_SHARED \
                     -I${QTINC}/mkspecs/linux-g++ -I${QTINC}/QtCore -I${QTINC}/QtNetwork  -I${QTINC} -c test_trackRT.cpp
	$(FTN) $(FFLAGS) test_trackRT.o -L/usr/local/lib -L${QTLIB} -lQtNetwork -lQtCore -lpthread  -lstdc++ $(SLIB) -o test_trackRT
	rm -f test_trackRT.o

# Make the Binary version (BNC 2.0/2.1 compatable version)
trackRTB : $(SLIB) trackRT.h trackRTObs.h
	g++ -c -pipe -O2 -Wall -W -D_REENTRANT -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_NETWORK_LIB -DQT_CORE_LIB -DQT_SHARED \
                     -I${QTINC}/mkspecs/linux-g++ -I${QTINC}/QtCore -I${QTINC}/QtNetwork  -I${QTINC} -c trackRTCommB.cpp
	$(FTN) $(FFLAGS) trackRTCommB.o -L/usr/local/lib -L${QTLIB} -lQtNetwork -lQtCore -lpthread  -lstdc++ $(SLIB) -o trackRTB
	rm -f trackRTCommB.o

$(TRACKRT_LIB) : \
	$(TRACKRT_LIB)(add_bfent.o) \
	$(TRACKRT_LIB)(earth_to_inert.o) \
	$(TRACKRT_LIB)(get_missrt.o) \
	$(TRACKRT_LIB)(obsrep.o) \
	$(TRACKRT_LIB)(remap_ambparn.o) \
	$(TRACKRT_LIB)(update_outfs.o) \
	$(TRACKRT_LIB)(assign_conrt.o) \
	$(TRACKRT_LIB)(form_ddrt.o) \
	$(TRACKRT_LIB)(get_sitert.o) \
	$(TRACKRT_LIB)(pred_step.o) \
	$(TRACKRT_LIB)(res_ambRT.o) \
	$(TRACKRT_LIB)(update_posrt.o) \
	$(TRACKRT_LIB)(atm_lib.o) \
	$(TRACKRT_LIB)(get_bflags.o) \
	$(TRACKRT_LIB)(init_estRT.o) \
	$(TRACKRT_LIB)(procobsRT.o) \
	$(TRACKRT_LIB)(theoryRT.o) \
	$(TRACKRT_LIB)(update_sp3.o) \
	$(TRACKRT_LIB)(comp_clockRT.o) \
	$(TRACKRT_LIB)(get_clk.o) \
	$(TRACKRT_LIB)(init_trackRT.o) \
	$(TRACKRT_LIB)(read_antmodrt.o) \
	$(TRACKRT_LIB)(TrackRTMods.o) \
	$(TRACKRT_LIB)(update_wlsstat.o) \
	$(TRACKRT_LIB)(comp_rawomc.o) \
	$(TRACKRT_LIB)(get_elev.o) \
	$(TRACKRT_LIB)(intantmod.o) \
	$(TRACKRT_LIB)(SaveObsA.o) \
	$(TRACKRT_LIB)(read_batrt.o)
	ranlib trackRT_lib.a

$(TRACKRT_LIB)(add_bfent.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(earth_to_inert.o): ../includes/const_param.h
$(TRACKRT_LIB)(get_missrt.o): ../includes/const_param.h trackRT_cmds.h trackRT.h
$(TRACKRT_LIB)(obsrep.o): obsInternal.h
$(TRACKRT_LIB)(remap_ambparn.o): ../includes/const_param.h trackRT.h
$(TRACKRT_LIB)(update_outfs.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(assign_conrt.o): trackRT_cmds.h
$(TRACKRT_LIB)(form_ddrt.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(get_sitert.o): ../includes/const_param.h trackRT_cmds.h trackRT.h
$(TRACKRT_LIB)(pred_step.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(res_ambRT.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(update_posrt.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(atm_lib.o): ../includes/const_param.h trackRT.h
$(TRACKRT_LIB)(get_bflags.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(init_estRT.o): trackRT.h trackRTObs.h
$(TRACKRT_LIB)(procobsRT.o): trackRT.h trackRTObs.h
$(TRACKRT_LIB)(theoryRT.o): ../includes/const_param.h trackRT.h
$(TRACKRT_LIB)(update_sp3.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(comp_clockRT.o): trackRT.h
$(TRACKRT_LIB)(get_clk.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(init_trackRT.o): ../includes/const_param.h trackRT.h
$(TRACKRT_LIB)(read_antmodrt.o): ../includes/const_param.h trackRT.h
$(TRACKRT_LIB)(TrackRTMods.o): trackRT.h trackRTObs.h
$(TRACKRT_LIB)(SaveObsA.o): trackRT.h trackRTObs.h
$(TRACKRT_LIB)(update_wlsstat.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(comp_rawomc.o): ../includes/const_param.h trackRT.h trackRTObs.h
$(TRACKRT_LIB)(get_elev.o): ../includes/const_param.h
$(TRACKRT_LIB)(read_batrt.o): trackRT_cmds.h trackRT.h

.c.a:
	$(CC) -c $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

# End of Makefile
