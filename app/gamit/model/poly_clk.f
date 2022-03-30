      subroutine poly_clk( ipass )
c
c PURPOSE: subroutine to compute receiver clock correction from
c          polynomial coefficients
c
c PARAMETERS:
c         IN: ipass  : clock iteration pass counter                      I*4
c           ** all these in model.h **
c             klock  : type of receiver clock correction 1=none
c                      2=polynomial, 3=polynomial/epoch                  I*4
c             jd0    : julian day of initial epoch and rec poly ref time I*4
c             t0     : time of day in secs (UTC) of initial epoch and
c                      receiver polynomial reference time                R*8
c             jdobs  : julian day of current epoch (no rec clock corrn)  I*4
c             tobs   : time of day in secs UTC of current epoch
c                      (no receiver clock correction applied)            R*8
c             clkepc : receiver clock polynomial offset coefficient      R*8
c             clkrat : receiver clock polynomial rate coefficient        R*8
c             clkacc : receiver clock polynomial acceleration coeff      R*8
c             clkcub : receiver clock polynomial cubic coefficient       R*8
c             rclock0: rclock derived from epoch wise approach from
c                      pass one of clock iteration                       R*8
c
c        OUT: rclock : receiver clock offset term from UTC time, used to
c                      correct jdobs and tobs to true observation time   R*8
c
c SUBROUTINES CALLED: timdif,report_stat
c
c CREATED: 27th DEC 1993               LAST MODIFIED: 27th DEC 1993
c
c AUTHOR: rwk, put in SR by S McClusky.
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1993
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/model.h'                


      integer*4 ipass
c
      real*8 tr,timdif

      if (klock .eq. 1) then
c        apply no correction at all
         rclock = 0.d0
         clkepc = 0.d0
         clkrat = 0.d0
         clkacc = 0.d0
      else if (klock .eq. 2) then
c        from input polynomial coefficients in s-file.
         tr = timdif(jdobs,tobs,jd0,t0)
         rclock = clkepc
     .          + clkrat*tr
     .          + clkacc*tr*tr
     .          + clkcub*tr*tr*tr
      else if (klock .eq. 3) then
c        use epoch by epoch correction, but use poly.
c        coefficients in S-file as first guess
         if (ipass .eq. 1) then
            tr = timdif(jdobs,tobs,jd0,t0)
            rclock = clkepc
     .             + clkrat*tr
     .             + clkacc*tr*tr
     .             + clkcub*tr*tr*tr
         else
            rclock = rclock0
         endif
      else
         call report_stat('FATAL','MODEL','poly_clk',klock,
     .   'Error, unknown receiver clock model type (klock): ',0)
      endif           
cd      print *,'POLY_CLK klock clkepc clkrat tr rclock0 rclock '
cd     .                ,  klock,clkepc,clkrat,tr,rclock0,rclock
      return
      end
