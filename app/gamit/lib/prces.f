Copyright (c) Massachusetts Institute of Technology and the Regents of the University of
California at San Diego, 1995. All rights reserved.

      Subroutine PRCES( fjd,oblq,prec,frame,precmod )

C     Form the precession matrix.  Currently implemented for 1950.0 from
c     either the IAU68 and IAU76 values of the precession formulas.  Prior
c     to May 1995, the IAU68 values were hard wired.  After that date,
c     the computations are controlled by the variables "frame" and "precmod"
c     to determine whether J2000 or B1950 inertial frame is to be used, and
c     whether IAU68 or IAU76 precession model is to be used
c
c     ARC can process in any of the above combinations, MODEL will use whatever
c     was used in creating the tfile, the orbit programs will use either whatever
c     is in the tfile (inert-efixed) or ask the user which frame is required for
c     the inertial frame (efixed-inertial). In this case ,precession models will
c     be assigned as follows: B1950 -> IAU68 ; J2000 -> IAU76
c
c     The old expressions are from the Explanatory Supplement to the American
c     Ephemeris and Nautical Almanac (AENA) 1968.  The new formulas are from
c     from Lieske et al., Astron. Astrophys. 58, 1-16, 1977.

c     The sense of the matrix is from the reference mean epoch to the
c     current epoch.

c     This version by R. King (Feb 92) and P. Tregoning (Apr 95) from earlier
c     versions by R. Abbot and Y. Bock.

c     Input:

c       FJD      Epoch of Date (PEP JD = true JD + 0.5)
c       frame    Inertial reference frame
c       precmod  Precession model to be used
c
c      Output:

c       PREC  Precession matrix
c       OBLQ  Mean obliquity of date

      implicit none

      character*5 frame, precmod
      character*80 prog_name

      real*8 fjd,oblq,prec(3,3),ttrm(3),z(3)
     .     , ocof(3),zcof(3,3),oblq0,oblq0x
     .     , casr,pi,ts,tc,s,t,dot
     .     , pzeta(3,3),pz(3,3),ptheta(3,3),work(3,3)

      integer*4 len,rcpar

      logical first

c     This used only for IAU76 expressions
      data  first/.true./

      save oblq0,ocof,zcof,first

c     get the calling module for report_stat
      len = rcpar(0,prog_name)

      pi=4.d0*datan(1.d0)
      casr=pi/180.d0/3600.D0


c  Compute T and t in precession formulas
c    T (=TC) is the time interval in Julian centuries of TDB between
c            2000.0 and 1950.0
c    t (=TS) is the time interval in Julian centuries of TDB between
c            epoch of date and reference epoch (B1950 or J2000)
c
c  If we want J2000 then just make TC = 0 (default J2000)
      tc = 0.d0
      ts = (fjd - 0.5d0 - 2451545.d0)/36525.d0
      if(frame.eq.'B1950')then
        tc = (2433282.423D0 - 2451545.D0)/36525.D0
        ts = (fjd - 0.5D0 - 2433282.423D0)/36525.d0
      elseif( frame.eq.'J2000') then
        tc = 0.d0
        ts = (fjd - 0.5d0 - 2451545.d0)/36525.d0
      else
        call report_stat('FATAL',prog_name,'lib/prces',frame
     .                  ,'Reference frame not recognized:',0)
      endif

c      print *,'PRCES frame precmod tc ',frame,precmod,tc

c  catch the case where we have J2000 and IAU68
      if(frame.eq.'J2000'.and.precmod.eq.'IAU68') then
         call report_stat('FATAL',prog_name,'lib/prces',' '
     .                   ,'J2000 and IAU68 are incompatible',0)
      endif

c      print *,'PRCES: frame and model ',frame,precmod
c  determine which precession model to use
      if(precmod.eq.'IAU68' ) then

c**     Code for the Old system (Ref: AENA, 1963)

        S = TS*36525.d0 + 1.8262423D+04
C       T is tropical centuries since fixed epoch 1950.0
        T = TS*36525.d0/36524.21988D0
        OBLQ = 4.093197552D-01 - S*(6.21795945D-09 + S*(2.1441068743D-17
     $     - S*1.800871677D-22))
        z(1) = t*(2.304948D+03 + t*( 0.302D0 + t*0.179D-01))*casr
        z(2) = t*(2.304948D+03 + t*( 1.093D0 + t*0.192D-01))*casr
        z(3) = t*(2.004255D+03 + t*(-0.426D0 - t*0.416D-01))*casr

      else

c ***   Code for the IAU76 system (Ref: Lieske et al., Astron. Astrophys. 58, 1-16, 1977)
        ttrm(1) = casr*ts
        ttrm(2) = ttrm(1)*ts
        ttrm(3) = ttrm(2)*ts
C
        if ( first) then
           oblq0 = 84381.448D0 - 46.8150D0*tc - 0.00059D0*tc**2
     .             + 0.001813D0*tc**3
C          At 1950.0 this give OBLQ0 = 84404.85522
           ocof(1) = -46.8150D0 - 0.00117D0*tc + 0.005439D0*tc**2
           ocof(2) = -0.00059D0 + 0.005439D0*tc
           ocof(3) =  0.001813D0
           zcof(1,1) = 2306.2181D0 + 1.39656D0*tc - 0.000139D0*tc**2
C          (UNSW Mon. 9 indicates -0.000344 (?))
           zcof(1,2) = 0.30188D0 - 0.000345D0*tc
           zcof(1,3) = 0.017998D0
           zcof(2,1) = zcof(1,1)
           zcof(2,2) = 1.09468D0 + 0.000066D0*tc
           zcof(2,3) = 0.018203D0
           zcof(3,1) = 2004.3109D0 - 0.85330D0*tc - 0.000217D0*tc**2
           zcof(3,2) = -0.42665D0 - 0.000217D0*tc
           zcof(3,3) = -0.041833D0
           first = .false.
        endif

        oblq0x=oblq0
        oblq=oblq0x*casr + dot(ocof,ttrm)
        call matmpy(zcof,ttrm,z,3,3,1)
C          zeta = z(1),  z = z(2),  theta = z(3)


      endif

c      print *,'PRCES: precmod tc ts zcof oblq ',precmod,tc,ts,zcof,oblq

c  Calculate the precession matrix from the reference to the current epoch

c     Reference:  Mueller, p. 65
      Call ROTMAT(-z(1),3,pzeta)
      Call ROTMAT(z(3),2,ptheta)
      Call ROTMAT(-z(2),3,pz)
      Call MATMPY(ptheta,pzeta,work,3,3,3)
      Call MATMPY(pz,work,prec,3,3,3)

      return
      end
