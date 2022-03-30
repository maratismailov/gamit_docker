c
      SUBROUTINE GEOXYZ(semi,finv,tx,ty,tz,ALAT,ALONG,hght,
     1    x,y,z,iflag,HEIGHT)
C
C     IFLAG=1: CONVERTS GEODETIC COORDINATES TO CARTESIAN COORDINATES
C     IFLAG=2: CONVERTS CARTESIAN CCORDINATES TO GEODETIC COORDINATES
c     alat alon unit: radian
c     semi, height, hght: meter
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer iflag
C
      if (semi.lt.6.0d6.or.finv.lt.2.5d2) then
         print*,' Wrong semi and finv! ',semi,finv
         stop
      endif
      TWOPI=8.D0*DATAN(1.D0)
      F=1.D0/finv
      E2=2.D0*F-F*F
      IF(IFLAG.EQ.2) GO TO 10
      SINLAT=DSIN(ALAT)
      COSLAT=DCOS(ALAT)
      SINLON=DSIN(ALONG)
      COSLON=DCOS(ALONG)
      CURVN=semi/(DSQRT(1.D0-E2*SINLAT*SINLAT))
C
      X=(CURVN+hght)*COSLAT*COSLON+tx
      Y=(CURVN+hght)*COSLAT*SINLON+ty
      Z=(CURVN*(1.D0-E2)+hght)*SINLAT+tz
C
      GO TO 20
C
   10 CONTINUE
      x = x-tx
      y = y-ty
      z = z-tz
      ALONG=DATAN2(Y,X)
      IF(ALONG.LT.0.D0) ALONG=ALONG+TWOPI
C     STARTING VALUE FOR LATITUDE ITERATION
      SQR=DSQRT(X*X+Y*Y)
      ALAT0=DATAN2(Z/SQR,1.D0-E2)
      ALAT=ALAT0
  40  SINLAT=DSIN(ALAT)
      CURVN=semi/(DSQRT(1.D0-E2*SINLAT*SINLAT))
      ALAT=DATAN2((Z+E2*CURVN*SINLAT),SQR)
C     ITERATE TO THE MILLIMETER LEVEL
      IF(DABS(ALAT-ALAT0).LT.1.D-10) GO TO 30
      ALAT0=ALAT
      GO TO 40
   30 CONTINUE
      CUTOFF=80.D0*TWOPI/360.D0
      IF(ALAT.GT.CUTOFF) GO TO 50
      hght=(SQR/DCOS(ALAT))-CURVN
      GO TO 20
   50 hght=Z/DSIN(ALAT)-CURVN+E2*CURVN
C
C  FOR JACOBIAN (GEODETIC TO XYZ)
 20   HEIGHT=CURVN+hght

      RETURN
      END
