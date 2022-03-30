Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.

      Subroutine SIDMAT( jd,fract,ut1utc,EqE,sidten,sdtmat,sidtm
     .                 , precmod )

c     Compute the Sidereal time rotation matrix
C

C       Parameters

c   Input:

C      JD        Julian day (PEP JD)
C      FRACT     UTC fraction of day
C      UT1UTC    UT1-UTC
C      EqE       Equation of Equinoxes (dpsi*cos(eps))
c      precmod - precession model

c   Output:

C      SIDTEN    Sidereal time rotation matrix
C      SDTMAT    Sidereal velocity matrix
c      SIDTM    Sidereal time (radians)

C
      implicit none

      character*5 precmod

      integer*4 jd

      real*8 fract,ut1utc,EqE,sidtm,sidten,sdtmat,sidtm0,sidvel
     .     , stheta,ctheta,twopi,ut1

      dimension sidten(3,3),sdtmat(3,3)

C
      twopi=8.d0*datan(1.d0)

      call sidtim( jd,sidtm0,sidvel,precmod )
c     units are radians, radians/s
                                  
c      write(*,100)sidtm0
c100   format('sidtm0: ',f22.14)

C UT1 in seconds of time

      ut1 = fract*86400.d0 + ut1utc


C GAST = GMST0 + sidvel*UT1 + Eq. E

      sidtm=sidtm0 + ut1*sidvel + EqE
C
      sidtm=dmod(sidtm,twopi)


C     CALCULATE SIDEREAL TIME MATRIX
        CTHETA = DCOS(SIDTM)
        STHETA = DSIN(SIDTM)
        SIDTEN(1,1) =  CTHETA
        SIDTEN(1,2) =  STHETA
        SIDTEN(1,3) =  0.D0
        SIDTEN(2,1) = -STHETA
        SIDTEN(2,2) =  CTHETA
        SIDTEN(2,3) =  0.D0
        SIDTEN(3,1) =  0.D0
        SIDTEN(3,2) =  0.D0
        SIDTEN(3,3) =  1.D0

C Compute sidereal velocity (in radians/sec)
C
        SDTMAT(1,1) = -SIDVEL*STHETA
        SDTMAT(1,2) = +SIDVEL*CTHETA
        SDTMAT(1,3) =  0.0D0
        SDTMAT(2,1) = -SIDVEL*CTHETA
        SDTMAT(2,2) = -SIDVEL*STHETA
        SDTMAT(2,3) = 0.0D0
        SDTMAT(3,1) = 0.0D0
        SDTMAT(3,2) = 0.0D0
        SDTMAT(3,3) = 0.0D0
C
        RETURN
        END
