      SUBROUTINE SIDTIM(JD,SIDTM0,SIDVEL,precmod )
C
C     Compute Sidereal Time given a Julian Day Number
C     M.E.Ash  JULY 1964
C
C      JD=      Input Julian Day number
C      SIDTM0 = Output mean sidereal time at 0 hr UTC
C             = Mean Sidereal Time at Julian Date JD - . 5 (in radians)
C      SIDVEL = Output ratio (Mean Sideral Time)/(Universal Time)
C                 multiplied by twopi divided by 86,400
c      precmod -  precession model
C
      implicit none

      character*5 precmod

      INTEGER*4 JD

      real*8 twopi,sidtm0,sidvel,t,tu2000


      TWOPI= 8.D0*DATAN(1.D0)
C
      T= DBLE(JD-2415020)
      T=T-0.5D0

c Define TU2000 to be number of julian days since J2000 (in centuries)
c Extra 0.5 days subtracted to align PEP and Julian Days
       tu2000 = (dble(jd) - 2451545.d0 - 0.5d0)/36525.d0
c        write(*,'(a,1x,a5)')' in sidtim precmod is',precmod

      if ( precmod.eq.'IAU68' ) then

C.....Old expressions for 'SIDTM0' AND 'SIDVEL'
c
c PT 950424: old expressions in seconds:
c
c  sidvel = 1.0027379092650 + 5.8899998d-11*T
c            h  m  s
c  sidtm0 = 6 38 45.83600 + 8640184.5419964*T + 2.5434635d-6*T**2
c           where T is time since 1900 (in centuries)
c
        SIDVEL =7.2921158546827D-5+T*1.1727115D-19
        SIDTM0=1.73993589372D0+T*(1.7202791266D-2+T*5.06409D-15)
c        print*,' using old sidtim '
      else
c        print*,' using new sidtim '
C.....The following equations for 'SIDTM0' and 'SIDVEL' should be
C.....used with the IAU76 precession
C.....The new defination was recommended by IAU Commisions 4, 19,
C.....and 31 so that there will be no change in the rate of UT1
C.....when the IAU(76) precession constant is used, or a change
C.....in the value of UT1 when the FK5 equinox is used.
c
c PT 950424: code in radians changed to IERS code in seconds of time
c        SIDVEL = 7.2921158548774D-5 + T*(1.1750609D-19 - T*3.22D-28)
c       write(*,2314)' gamit ',sidvel*86400.d0/twopi
c2314   format(a7,f20.18)
c        SIDTM0=1.739935359516D0+T*(1.72027914345D-2+T*(5.076246D-15
c     $     - T*9.25D-24))
c
c  code in seconds was:
c
c  sidvel = 1.002 737 909 291 8 + 5.90180d-11*T - 1.617260d-19*T**3
c
c            h  m  s                   s
c  sidtm0 = 6 38 45.828654203 + 8640184.6266264*T + 2.5495689357841D-06*T**2 - 4.6458569297080D-15*T**3
c           where T is time since 1900 (in centuries)
c
c IERS formula is: (ref. IERS technical note #13, Aoki et al. Astron, Astrophys, 105, 359-361, 1982.)
c
c  sidvel = 1.002 737 909 350 795 + 5.9006d-11*Tu - 5.9d-15*Tu**2
c            h  m  s               s             s                s
c  sidtm0 = 6 41 50.54841 + 8640184.812866*Tu + 0.093104*Tu**2 - 6.2d-6*Tu**3
c           where Tu is time since J2000 (in centuries)
c
c day 076,1995 difference between formulations
c
c SIDTM0:  1.2036D-5 seconds
c SIDVEL:  2.6880053d-11
c
c  USE the IERS coding

c  PT 950424: Code the iers standard for sidtm0 in seconds of time.
        sidvel=1.002737909350795d0 +5.9006d-11*tu2000 -5.9d-15*tu2000**2
c       write(*,2314),' IERS  ',sidvel
        sidtm0 =  24110.54841d0
     .       + 8640184.812866d0*tu2000 + 0.093104d0*tu2000**2
     .        - 6.2d-6*tu2000**3

c  convert seconds to radians
        sidtm0 = sidtm0*twopi/86400.d0
        sidvel = sidvel*twopi/86400.d0

      endif
C
      SIDTM0=DMOD(SIDTM0,TWOPI)
c      write(*,900)precmod,t,tu2000,sidtm0,sidvel
c900   format('In SIDTIM: ',a5,1x,f26.18,1x,f22.18,1x,f22.18,1x,f22.18)
C
c      print*,' in sidtim sidtm0 is ',sidtm0
c      stop
      RETURN
      END
