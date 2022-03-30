# Makefile for plot on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:41:23 UTC 2021

all: plot_lib cplotx multiplot multibase



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

PLOT_LIB = plot_lib.a
SLIB = ./plot_lib.a x1000/x1000.a ../gen_util/gen_util_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

plot_lib : $(PLOT_LIB)

cplotx : cplot.f  $(SLIB) plot_bd.f plot_com.h plot_ema.h plot_param.h
	$(FTN) $(FFLAGS) cplot.f  $(SLIB) -L$(X11LIB) -lX11 -o cplotx
	rm -f cplot.o

multiplot : multiplot.f  $(SLIB) 
	$(FTN) $(FFLAGS) multiplot.f  $(SLIB) -L$(X11LIB) -lX11 -o multiplot
	rm -f multiplot.o

multibase : multibase.f  $(SLIB) 
	$(FTN) $(FFLAGS) multibase.f  $(SLIB) -L$(X11LIB) -lX11 -o multibase
	rm -f multibase.o

$(PLOT_LIB) : \
	$(PLOT_LIB)(plot_spec.o) \
	$(PLOT_LIB)(plts1.o) \
	$(PLOT_LIB)(plts2.o) \
	$(PLOT_LIB)(plts3.o) \
	$(PLOT_LIB)(plts8.o) \
	$(PLOT_LIB)(plts9.o) \
	$(PLOT_LIB)(pltsa.o) \
	$(PLOT_LIB)(pltsp.o) \
	$(PLOT_LIB)(setup.o) \
	$(PLOT_LIB)(finish_up.o) \
	$(PLOT_LIB)(report.o) \
	$(PLOT_LIB)(command_erro.o) \
	$(PLOT_LIB)(get_cmd.o) \
	$(PLOT_LIB)(trim_lead.o) \
	$(PLOT_LIB)(get_int.o) \
	$(PLOT_LIB)(get_xy_scale.o) \
	$(PLOT_LIB)(draw.o) \
	$(PLOT_LIB)(draw_errb.o) \
	$(PLOT_LIB)(get_file.o) \
	$(PLOT_LIB)(read_data.o) \
	$(PLOT_LIB)(check_bound.o) \
	$(PLOT_LIB)(get_view.o) \
	$(PLOT_LIB)(set_view.o) \
	$(PLOT_LIB)(get_scale.o) \
	$(PLOT_LIB)(set_scale.o) \
	$(PLOT_LIB)(dxaxis_mn.o) \
	$(PLOT_LIB)(format_label.o) \
	$(PLOT_LIB)(label_ax1.o) \
	$(PLOT_LIB)(label_ax0.o) \
	$(PLOT_LIB)(draw_axis.o) \
	$(PLOT_LIB)(axis_cmd.o) \
	$(PLOT_LIB)(com_space.o) \
	$(PLOT_LIB)(time_space.o) \
	$(PLOT_LIB)(first_tic1.o) \
	$(PLOT_LIB)(first_tic0.o) \
	$(PLOT_LIB)(inc_cax.o) \
	$(PLOT_LIB)(draw_tic.o) \
	$(PLOT_LIB)(axis_setup.o) \
	$(PLOT_LIB)(conv_mmtow.o) \
	$(PLOT_LIB)(dxaxis_mx.o) \
	$(PLOT_LIB)(dyaxis_mn.o) \
	$(PLOT_LIB)(dyaxis_mx.o) \
	$(PLOT_LIB)(label_axis.o) \
	$(PLOT_LIB)(put_label.o) \
	$(PLOT_LIB)(read_label.o) \
	$(PLOT_LIB)(set_charsz.o) \
	$(PLOT_LIB)(report_error.o) \
	$(PLOT_LIB)(status.o) \
	$(PLOT_LIB)(report_scale.o) \
	$(PLOT_LIB)(out_names.o) \
	$(PLOT_LIB)(get_charsz.o) \
	$(PLOT_LIB)(write_label.o) \
	$(PLOT_LIB)(report_windo.o) \
	$(PLOT_LIB)(set_texl.o) \
	$(PLOT_LIB)(process.o) \
	$(PLOT_LIB)(s2mov.o) \
	$(PLOT_LIB)(sr2mv.o) \
	$(PLOT_LIB)(scj.o) \
	$(PLOT_LIB)(save_label.o) \
	$(PLOT_LIB)(output_label.o) \
	$(PLOT_LIB)(swh.o) \
	$(PLOT_LIB)(rotor.o) \
	$(PLOT_LIB)(comlf.o) \
	$(PLOT_LIB)(draw_map.o) \
	$(PLOT_LIB)(get_poly_win.o) \
	$(PLOT_LIB)(poly_fit.o) \
	$(PLOT_LIB)(pdraw.o) \
	$(PLOT_LIB)(mark_window.o) \
	$(PLOT_LIB)(identify.o) \
	$(PLOT_LIB)(make_line.o) \
	$(PLOT_LIB)(err_lib.o) \
	$(PLOT_LIB)(advance_rc.o) \
	$(PLOT_LIB)(network.o) \
	$(PLOT_LIB)(strain_lib.o)
	ranlib plot_lib.a

$(PLOT_LIB)(plot_spec.o): plot_com.h plot_param.h
$(PLOT_LIB)(plts1.o): plot_com.h plot_ema.h plot_param.h
$(PLOT_LIB)(plts2.o): plot_com.h plot_ema.h plot_param.h
$(PLOT_LIB)(plts3.o): plot_com.h plot_param.h
$(PLOT_LIB)(plts8.o): plot_com.h plot_param.h
$(PLOT_LIB)(plts9.o): plot_com.h plot_ema.h plot_param.h
$(PLOT_LIB)(pltsa.o): plot_com.h plot_param.h
$(PLOT_LIB)(pltsp.o): plot_com.h plot_ema.h plot_param.h
$(PLOT_LIB)(setup.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_xy_scale.o): plot_com.h plot_param.h
$(PLOT_LIB)(draw.o): plot_com.h plot_param.h
$(PLOT_LIB)(draw_errb.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_file.o): plot_com.h plot_param.h
$(PLOT_LIB)(read_data.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_view.o): plot_com.h plot_param.h
$(PLOT_LIB)(set_view.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_scale.o): plot_com.h plot_param.h
$(PLOT_LIB)(set_scale.o): plot_com.h plot_param.h
$(PLOT_LIB)(dxaxis_mn.o): plot_com.h plot_param.h
$(PLOT_LIB)(draw_axis.o): plot_com.h plot_param.h
$(PLOT_LIB)(time_space.o): plot_com.h plot_param.h
$(PLOT_LIB)(draw_tic.o): plot_com.h plot_param.h
$(PLOT_LIB)(dxaxis_mx.o): plot_com.h plot_param.h
$(PLOT_LIB)(dyaxis_mn.o): plot_com.h plot_param.h
$(PLOT_LIB)(dyaxis_mx.o): plot_com.h plot_param.h
$(PLOT_LIB)(label_axis.o): plot_com.h plot_param.h
$(PLOT_LIB)(put_label.o): plot_com.h plot_param.h
$(PLOT_LIB)(read_label.o): plot_com.h plot_param.h
$(PLOT_LIB)(set_charsz.o): plot_com.h plot_param.h
$(PLOT_LIB)(status.o): plot_com.h plot_param.h
$(PLOT_LIB)(report_scale.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_charsz.o): plot_com.h plot_param.h
$(PLOT_LIB)(write_label.o): plot_com.h plot_param.h
$(PLOT_LIB)(report_windo.o): plot_com.h plot_param.h
$(PLOT_LIB)(set_texl.o): plot_com.h plot_param.h
$(PLOT_LIB)(process.o): plot_com.h plot_param.h
$(PLOT_LIB)(s2mov.o): plot_com.h plot_param.h
$(PLOT_LIB)(sr2mv.o): plot_com.h plot_param.h
$(PLOT_LIB)(scj.o): plot_com.h plot_param.h
$(PLOT_LIB)(save_label.o): plot_com.h plot_param.h
$(PLOT_LIB)(output_label.o): plot_com.h plot_param.h
$(PLOT_LIB)(swh.o): plot_com.h plot_param.h
$(PLOT_LIB)(rotor.o): plot_com.h plot_param.h
$(PLOT_LIB)(comlf.o): plot_com.h plot_param.h
$(PLOT_LIB)(draw_map.o): plot_com.h plot_param.h
$(PLOT_LIB)(get_poly_win.o): plot_com.h plot_param.h
$(PLOT_LIB)(poly_fit.o): plot_com.h plot_param.h
$(PLOT_LIB)(pdraw.o): plot_com.h plot_param.h
$(PLOT_LIB)(mark_window.o): plot_com.h plot_param.h
$(PLOT_LIB)(identify.o): plot_com.h plot_param.h
$(PLOT_LIB)(err_lib.o): ../includes/const_param.h plot_com.h plot_param.h
$(PLOT_LIB)(network.o): plot_com.h plot_param.h
$(PLOT_LIB)(strain_lib.o): plot_com.h plot_param.h

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
