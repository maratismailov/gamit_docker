      real*8 function cfa(p,t,wetvar,h2otyp,phi,h,elev)
c
c     Calculates the CfA-2.2 mapping function (Davis et al., Radio
c     Science 20, 1593-1607, 1985)
c
c     INPUT:    p      Pressure, mbar
c               t      Temperature, deg C
c               wetvar Water vapor variable, defined below
c               h2otyp H2OTYP = 'R' --> WETVAR is relative humidity (0-1)
c                      H2OTYP = 'D' --> WETVAR is dew point temperature (deg C)
c               phi    Geocentric latitude (radians)
c               h      Elevation above ellipsoid, km
c               elev   Elevation angle, radians
c
      implicit none
      real*8 p, t, wetvar, phi, h, elev, e, rh
      real*8 a0, a1, a2, a3, a5
      real*8 b0, b1, b2, b3, b5
      real*8 c0
      real*8 a, b, c
      real*8 sinel, tanel
c
      character*1 h2otyp,UPPERC
c
      data a0 /  0.001185D+00/
      data a1 /  0.6071D-04 /
      data a2 / -0.1471D-03 /
      data a3 /  0.3072D-02 /
c      data a4 /  0.1965D-01 /
      data a5 / -0.5645D-02 /
c
      data b0 /  0.001144D+00 /
      data b1 /  0.1164D-04 /
      data b2 /  0.2795D-03 /
      data b3 /  0.3109D-02 /
c      data b4 /  0.3038D-01 /
      data b5 / -0.1217D-01 /
c
      data c0 / -0.009D+00 /
c
c.... Get partial pressure of water vapor
      if (h2otyp .eq. UPPERC('R')) then
        call wpress(1,wetvar,e,t)
c**        e = wpress(wetvar,t)
      elseif( htotyp .eq.'UPPER('D')) then
        call wpress(1,1.d0,wetvar,t)
c**        e = wpress(1.0D0,wetvar)
      end if
c
c.... The following expressions assume that the temperature lapse rate
c     is its nominal value of -6.5 K/km.  The tropopause height is assumed
c     to be the nominal sea-level value of 11.231 km less the height
c     of the station.
c

c** temporary fix for elev=90. case

      if( elev.eq.90.d0 ) then
        cfa = 1.d0
        goto 999
      endif

c.... Calculate sine and tangent of elevation
      sinel = sin(elev)
c** temporary fix for elev=90. case
      if( dabs(sinel-1.d0).lt.1.d-10 ) then
        cfa = 1.d0
        goto 999
      endif
      tanel = sinel / cos(elev)
c
c.... Caclulate the a, b, and c functions
      a = a0 * (1.0D+00 + a1 * (p -1000.0)
     .                  + a2 * e
     .                  + a3 * (t - 20.0)
     .                  + a5 * (-h)        )
c
      b = b0 * (1.0D+00 + b1 * (p - 1000.0)
     .                  + b2 * e
     .                  + b3 * (t - 20.0)
     .                  + b5 * (-h)        )
c
      c = c0
c
c.... Calculate mapping function
      cfa = 1.0D0 / (sinel + a / (tanel + b / (sinel + c)))
c
 999  return
      end
