Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved. 

      Subroutine GEOXYZ( semi,finv,alat,along,hght,geodrad,x,y,z,iflag )
C
C     iflag=1: Converts geodetic coordinates to Cartesian coordinates
c     iflag=2: Converst Cartesian coordinate to Geodetic coordinates

c        semi :  semi-major axis of reference ellipsoid
c        finv :  inverse flattening (1/f)
c        alat :  geodetic latitude in radians
c        along:  longitude in radians   
c        hght :  geodetic height
c        x,y,z:  Cartesian coordinates  
c        geodrad: geodetic height + inverse radius of (ellipsoid) curvature ('N' in Mueller),
c                 used to compute the Jacobian


      IMPLICIT none

      integer iflag
      real*8 semi,finv,alat,along,hght,x,y,z,geodrad,twopi,f,e2
     .       ,sinlat,sinlon,coslat,coslon,curvn,sqr,alat0,cutoff

      TWOPI=8.D0*DATAN(1.D0)
      F=1.D0/FINV
      E2=2.D0*F-F*F
      IF(IFLAG.EQ.2) GO TO 10
      SINLAT=DSIN(ALAT)
      COSLAT=DCOS(ALAT)
      SINLON=DSIN(ALONG)
      COSLON=DCOS(ALONG)
      CURVN=SEMI/(DSQRT(1.D0-E2*SINLAT*SINLAT))
C
      X=(CURVN+HGHT)*COSLAT*COSLON
      Y=(CURVN+HGHT)*COSLAT*SINLON
      Z=(CURVN*(1.D0-E2)+HGHT)*SINLAT
C
C  FOR JACOBIAN (GEODETIC TO XYZ)
      geodrad=CURVN+HGHT
C
      GO TO 20
C
   10 CONTINUE
      ALONG=DATAN2(Y,X)
      IF(ALONG.LT.0.D0) ALONG=ALONG+TWOPI
C     STARTING VALUE FOR LATITUDE ITERATION
      SQR=DSQRT(X*X+Y*Y)
C Changed by D.D.Dong:  Old:  ALAT0=DATAN2(1.D0-E2),Z/SQR)
      ALAT0=DATAN2(Z/SQR,1.D0-E2)
      ALAT=ALAT0
  40  SINLAT=DSIN(ALAT)
      CURVN=SEMI/(DSQRT(1.D0-E2*SINLAT*SINLAT)) 
      ALAT=DATAN2((Z+E2*CURVN*SINLAT),SQR)
C     ITERATE TO THE MILLIMETER LEVEL
      IF(DABS(ALAT-ALAT0).LT.1.D-10) GO TO 30
      ALAT0=ALAT
      GO TO 40
   30 CONTINUE
      CUTOFF=80.D0*TWOPI/360.D0
      IF(ALAT.GT.CUTOFF) GO TO 50
      HGHT=(SQR/DCOS(ALAT))-CURVN
      GO TO 20
   50 HGHT=Z/DSIN(ALAT)-CURVN+E2*CURVN 
      geodrad = 0.d0
   20 CONTINUE   
      RETURN
      END
