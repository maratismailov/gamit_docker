# Makefile for solve on Linux 5.10.83-1-lts by unimake at Fri Dec  3 17:24:04 UTC 2021

all: solve

FFLAGS = -O3 -Wuninitialized -Wunused -fno-f2c -ffast-math -m64

 

# ASSIGMENTS

# ASSIGMENTS
### LF95 flags ###
#FTN = lf95
#CC = gcc
#CFLAGS = -g
### IFC 7.0 flags ###
#FTN = ifc
#CC = icc
#CFLAGS = -O
### IFORT 8.0 flags ###
#FTN = ifort
#CC = icc
#CFLAGS = -O
### IFORT 9.0 flags ###
#FTN = ifort
#CC = icc
#CFLAGS = -O
### G77 flags ###
#FTN = g77
#CC = gcc
#CFLAGS = -g
### GCC V4 gfortran flags ###
FTN = gfortran
CC = gcc
#CFLAGS = -g
#### For 64-bit machines:
CFLAGS = -g -m64 -DADD64BIT

SOLVE_LIB = solve_lib.a
SLIB = ../lib/gamit_lib.a ../../libraries/comlib/com_lib.a 
X11LIB = /usr/lib/x86_64-linux-gnu
II = /usr/include

solve : solve.f $(SOLVE_LIB) $(SLIB) ../includes/dimpar.h models.h parameters.h solve.h
	$(FTN) $(FFLAGS) solve.f $(SOLVE_LIB) $(SLIB) -o solve
	rm -f solve.o

$(SOLVE_LIB) : \
	$(SOLVE_LIB)(abnew.o) \
	$(SOLVE_LIB)(add1.o) \
	$(SOLVE_LIB)(add2.o) \
	$(SOLVE_LIB)(add2i.o) \
	$(SOLVE_LIB)(addwgt.o) \
	$(SOLVE_LIB)(addwgtb.o) \
	$(SOLVE_LIB)(apply.o) \
	$(SOLVE_LIB)(atmcon.o) \
	$(SOLVE_LIB)(atmdif.o) \
	$(SOLVE_LIB)(atm_grad_part.o) \
	$(SOLVE_LIB)(atmstr.o) \
	$(SOLVE_LIB)(atpa.o) \
	$(SOLVE_LIB)(atpdal.o) \
	$(SOLVE_LIB)(avgload.o) \
	$(SOLVE_LIB)(avgvar.o) \
	$(SOLVE_LIB)(baslin.o) \
	$(SOLVE_LIB)(bcheck.o) \
	$(SOLVE_LIB)(bdeci.o) \
	$(SOLVE_LIB)(bexpec.o) \
	$(SOLVE_LIB)(bfwork.o) \
	$(SOLVE_LIB)(biset.o) \
	$(SOLVE_LIB)(bisopt.o) \
	$(SOLVE_LIB)(bnew.o) \
	$(SOLVE_LIB)(bnew1.o) \
	$(SOLVE_LIB)(bsort.o) \
	$(SOLVE_LIB)(bstat.o) \
	$(SOLVE_LIB)(check_adjust.o) \
	$(SOLVE_LIB)(copy1d.o) \
	$(SOLVE_LIB)(copy1i.o) \
	$(SOLVE_LIB)(covbln.o) \
	$(SOLVE_LIB)(covst1.o) \
	$(SOLVE_LIB)(covst2.o) \
	$(SOLVE_LIB)(cpmat2.o) \
	$(SOLVE_LIB)(csort.o) \
	$(SOLVE_LIB)(dbadst.o) \
	$(SOLVE_LIB)(dbar.o) \
	$(SOLVE_LIB)(dbcnt.o) \
	$(SOLVE_LIB)(dbias.o) \
	$(SOLVE_LIB)(dcheck.o) \
	$(SOLVE_LIB)(degfrm.o) \
	$(SOLVE_LIB)(degred.o) \
	$(SOLVE_LIB)(delcon.o) \
	$(SOLVE_LIB)(distnc.o) \
	$(SOLVE_LIB)(dopt.o) \
	$(SOLVE_LIB)(dot.o) \
	$(SOLVE_LIB)(double.o) \
	$(SOLVE_LIB)(dxyz_to_neu.o) \
	$(SOLVE_LIB)(eptime.o) \
	$(SOLVE_LIB)(filld.o) \
	$(SOLVE_LIB)(fillj.o) \
	$(SOLVE_LIB)(elem.o) \
	$(SOLVE_LIB)(filobs.o) \
	$(SOLVE_LIB)(filomc.o) \
	$(SOLVE_LIB)(filpar.o) \
	$(SOLVE_LIB)(fnddbi.o) \
	$(SOLVE_LIB)(formn1.o) \
	$(SOLVE_LIB)(formn2.o) \
	$(SOLVE_LIB)(fxbias.o) \
	$(SOLVE_LIB)(getcda.o) \
	$(SOLVE_LIB)(gethed.o) \
	$(SOLVE_LIB)(getslt.o) \
	$(SOLVE_LIB)(get_bias_scale.o) \
	$(SOLVE_LIB)(get_eop_apr.o) \
	$(SOLVE_LIB)(get_err_apr.o) \
	$(SOLVE_LIB)(get_grad_apr.o) \
	$(SOLVE_LIB)(get_narrowlane.o) \
	$(SOLVE_LIB)(get_sat_apr.o) \
	$(SOLVE_LIB)(get_sit_apr.o) \
	$(SOLVE_LIB)(get_solve_cut.o) \
	$(SOLVE_LIB)(get_zen_apr.o) \
	$(SOLVE_LIB)(get_widelane.o) \
	$(SOLVE_LIB)(glohed.o) \
	$(SOLVE_LIB)(gpgt.o) \
	$(SOLVE_LIB)(gradout.o) \
	$(SOLVE_LIB)(gtpgx.o) \
	$(SOLVE_LIB)(kepxyz.o) \
	$(SOLVE_LIB)(keywrd.o) \
	$(SOLVE_LIB)(lcloos.o) \
	$(SOLVE_LIB)(lcnorm.o) \
	$(SOLVE_LIB)(live.o) \
	$(SOLVE_LIB)(loos12.o) \
	$(SOLVE_LIB)(lrgint.o) \
	$(SOLVE_LIB)(lsqdo1.o) \
	$(SOLVE_LIB)(lsqerr.o) \
	$(SOLVE_LIB)(lsqint.o) \
	$(SOLVE_LIB)(lsqio.o) \
	$(SOLVE_LIB)(lsqprt.o) \
	$(SOLVE_LIB)(lsquar.o) \
	$(SOLVE_LIB)(lweop.o) \
	$(SOLVE_LIB)(lc_solution.o) \
	$(SOLVE_LIB)(lc_soln12.o) \
	$(SOLVE_LIB)(lwgrad.o) \
	$(SOLVE_LIB)(lwsat.o) \
	$(SOLVE_LIB)(lwstat.o) \
	$(SOLVE_LIB)(lwzen.o) \
	$(SOLVE_LIB)(matprt.o) \
	$(SOLVE_LIB)(nbias.o) \
	$(SOLVE_LIB)(nbias1.o) \
	$(SOLVE_LIB)(nbiasp.o) \
	$(SOLVE_LIB)(nbiasr.o) \
	$(SOLVE_LIB)(newfil.o) \
	$(SOLVE_LIB)(nm2131.o) \
	$(SOLVE_LIB)(nm2233.o) \
	$(SOLVE_LIB)(nm32.o) \
	$(SOLVE_LIB)(nmu11.o) \
	$(SOLVE_LIB)(normd.o) \
	$(SOLVE_LIB)(nrmsc2.o) \
	$(SOLVE_LIB)(nrmscl.o) \
	$(SOLVE_LIB)(nu23.o) \
	$(SOLVE_LIB)(oline.o) \
	$(SOLVE_LIB)(opera.o) \
	$(SOLVE_LIB)(opera2.o) \
	$(SOLVE_LIB)(pclock.o) \
	$(SOLVE_LIB)(pgrad.o) \
	$(SOLVE_LIB)(pseuwl.o) \
	$(SOLVE_LIB)(pzenth.o) \
	$(SOLVE_LIB)(qhead1.o) \
	$(SOLVE_LIB)(qhead2.o) \
	$(SOLVE_LIB)(qhead3.o) \
	$(SOLVE_LIB)(qhead4.o) \
	$(SOLVE_LIB)(raddeg.o) \
	$(SOLVE_LIB)(remenu.o) \
	$(SOLVE_LIB)(read_bfl.o) \
	$(SOLVE_LIB)(read_bfl_sess.o) \
	$(SOLVE_LIB)(read_biases.o) \
	$(SOLVE_LIB)(reodrd.o) \
	$(SOLVE_LIB)(round_bias.o) \
	$(SOLVE_LIB)(resolve_wl.o) \
	$(SOLVE_LIB)(rmjunk.o) \
	$(SOLVE_LIB)(set_para.o) \
	$(SOLVE_LIB)(solve1.o) \
	$(SOLVE_LIB)(solve2.o) \
	$(SOLVE_LIB)(solve3.o) \
	$(SOLVE_LIB)(solvlc.o) \
	$(SOLVE_LIB)(sort1i.o) \
	$(SOLVE_LIB)(sort4i.o) \
	$(SOLVE_LIB)(sortbl.o) \
	$(SOLVE_LIB)(sphxyz.o) \
	$(SOLVE_LIB)(stack.o) \
	$(SOLVE_LIB)(sversn.o) \
	$(SOLVE_LIB)(tapdrv.o) \
	$(SOLVE_LIB)(trans.o) \
	$(SOLVE_LIB)(transx.o) \
	$(SOLVE_LIB)(updat3.o) \
	$(SOLVE_LIB)(update.o) \
	$(SOLVE_LIB)(upi.o) \
	$(SOLVE_LIB)(upl.o) \
	$(SOLVE_LIB)(uporb.o) \
	$(SOLVE_LIB)(uselss.o) \
	$(SOLVE_LIB)(wbias.o) \
	$(SOLVE_LIB)(weop.o) \
	$(SOLVE_LIB)(wgrad.o) \
	$(SOLVE_LIB)(wl_fixed.o) \
	$(SOLVE_LIB)(wlmode.o) \
	$(SOLVE_LIB)(write_soln.o) \
	$(SOLVE_LIB)(wsat.o) \
	$(SOLVE_LIB)(wstat.o) \
	$(SOLVE_LIB)(wzen.o) \
	$(SOLVE_LIB)(xyzkep.o) \
	$(SOLVE_LIB)(zenout.o) \
	$(SOLVE_LIB)(zero1d.o) \
	$(SOLVE_LIB)(zero1i.o) \
	$(SOLVE_LIB)(zero2d.o) \
	$(SOLVE_LIB)(zero2i.o)
	ranlib solve_lib.a

$(SOLVE_LIB)(abnew.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(addwgt.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(addwgtb.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(apply.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(atmcon.o): ../includes/dimpar.h
$(SOLVE_LIB)(atmdif.o): ../includes/dimpar.h
$(SOLVE_LIB)(avgload.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(bcheck.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(bfwork.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(biset.o): ../includes/dimpar.h
$(SOLVE_LIB)(bisopt.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(bnew.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(bnew1.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(bstat.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(check_adjust.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(covbln.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(covst1.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(covst2.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(cpmat2.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(dbadst.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(dbar.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(dbcnt.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(dbias.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(dcheck.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(delcon.o): ../includes/dimpar.h
$(SOLVE_LIB)(dopt.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(double.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(eptime.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(filld.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(filobs.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(filomc.o): ../includes/dimpar.h models.h solve.h
$(SOLVE_LIB)(filpar.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(fnddbi.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(formn1.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(formn2.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(fxbias.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(getcda.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(gethed.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(getslt.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_bias_scale.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_eop_apr.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(get_err_apr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_grad_apr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_narrowlane.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_sat_apr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_sit_apr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_solve_cut.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_zen_apr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(get_widelane.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(glohed.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(gradout.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(keywrd.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lcloos.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lcnorm.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(live.o): ../includes/dimpar.h
$(SOLVE_LIB)(loos12.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lrgint.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(lsqdo1.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lsqerr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lsqint.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lsqio.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lsqprt.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(lsquar.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lweop.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lc_solution.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lc_soln12.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lwgrad.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lwsat.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lwstat.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(lwzen.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(nbias.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(nbias1.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(nbiasp.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(nbiasr.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(nm2131.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(nm2233.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(nm32.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(nmu11.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(normd.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(nrmscl.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(nu23.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(oline.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(opera.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(opera2.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(pclock.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(pgrad.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(pseuwl.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(pzenth.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(qhead1.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(qhead2.o): ../includes/dimpar.h models.h solve.h
$(SOLVE_LIB)(qhead3.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(qhead4.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(remenu.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(read_bfl.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(read_bfl_sess.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(read_biases.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(reodrd.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(round_bias.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(resolve_wl.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(rmjunk.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(set_para.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(solve1.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(solve2.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(solve3.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(solvlc.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(transx.o): ../includes/dimpar.h
$(SOLVE_LIB)(updat3.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(update.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(upi.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(upl.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(uporb.o): ../includes/dimpar.h models.h parameters.h solve.h
$(SOLVE_LIB)(uselss.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(wbias.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(weop.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(wgrad.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(wlmode.o): ../includes/dimpar.h solve.h
$(SOLVE_LIB)(write_soln.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(wsat.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(wstat.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(wzen.o): ../includes/dimpar.h parameters.h solve.h
$(SOLVE_LIB)(zenout.o): ../includes/dimpar.h parameters.h solve.h

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
