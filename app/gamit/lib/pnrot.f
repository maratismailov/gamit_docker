Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.

      SUBROUTINE PNROT(INUT,JD,T,TDTGPST,EQE,PREC,RNUT,frame
     .                 ,precmod)
C
C Written by Yehuda Bock (1987)
C
C Compute the precession and nutation matrices

c     Input:
c       inut       Unit number for nutation table
c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
c       t          Seconds-of-day  (GPST)
c       tdtgpst    TDT - GPST  (seconds)
c       frame      Inertial frame
c       precmod    precession model
c
c      Output:
c        EqE          Equation of the equinoxes (rad)
c        prec         Precession matrix (rad)
c        rnut         Nutation matrix (rad)

      implicit none

      character*5 frame, precmod

      integer*4 jd,inut

      real*8 t,tdtgpst,eqe,prec,rnut,deps,oblq,dpsi,fjdtdt
C
      dimension prec(3,3),rnut(3,3)

c Convert GPS time to Terrestrial Dynamical Time (still PEP JD = true JD + 0.5)

      FJDTDT= DBLE(JD) + T/86400.D0 + TDTGPST/86400.d0

C Form the precession matrix


c     temporary to get 1950 to 2000 matrix
c      fjdsave = fjdtdt
c      fjdtdt = 2433282.923d0
c      call PRCES(FJDTDT,OBLQ,PREC,frame,precmod)
c      print *,'B1950 PREC: ',prec
c      fjdtdt = fjdsave

      call PRCES(FJDTDT,OBLQ,PREC,frame,precmod)

C Form the nutation matrix

      Call NUTTAB(inut,FJDTDT,OBLQ,DPSI,DEPS,RNUT)

C Compute the equation of the equinoxes for computation of GAST

      EqE=dpsi*dcos(oblq)

      return
      end
