      real*8 function marini(p,t,wetvar,h2otyp,phi,h,elev)
c
c     Calculates the Marini mapping function.
c
c     INPUT:
c       P         Total pressure, mbars
c       T         Temperature, deg C
c       WETVAR    Water vapor variable, defined by H2OTYP
c       H2OTYP    Defines WETVAR.  H2OTYP = 'R' indicates that WETVAR is relati
c                 humidity (0-1).  H2OTYP = 'D' indicates that WETVAR is the de
c                 point temperature (deg C).
c       PHI       Geocentric latitude, radians
c       H         Elevation above the geoid, km
c       ELEV      Pointing elevation, radians
c
      real*8 p, t, wetvar, phi, h, elev, beta, sinel
      real*8 A, B, ffun, saaszd, beta1

      character*1 h2otyp
c
c.... Calculate A and B values
      A = saaszd(p,t,wetvar,h2otyp,phi,h)
      B = 2.644D-03 * exp(-0.14372D+00 * h) / ffun(phi,h)
c
c.... Mapping function depends only on ratio
      beta  = B / A
      beta1 = 1.0D+00 + beta
c
c.... Sine of elevation
      sinel = sin(elev)
c
c.... Calculate mapping function
      marini = beta1 / (sinel + beta / beta1 / (sinel + 0.15D-01))
c
      end
