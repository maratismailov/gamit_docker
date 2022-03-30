C
      SUBROUTINE CARSPH(X,ALAT,ALONG,RAD)
C
C  Converts cartesian x,y,z to spherical coordinates
C
      implicit none
      real*8 x,twopi,rad,along,alat
      DIMENSION X(3)
C
      TWOPI= 8.D0*DATAN(1.D0)
C
      RAD= DSQRT(X(1)**2+X(2)**2+X(3)**2)
      ALONG= DATAN2(X(2),X(1))
C      IF( ALONG.LT.0.D0 ) ALONG=ALONG+TWOPI
      ALAT= DASIN(X(3)/RAD)
      ALONG= ALONG*360.D0/TWOPI
      ALAT=  ALAT*360.D0/TWOPI
C
      RETURN
      END
