      SUBROUTINE HELINV(GLAT1,E1,GLAT2,E2,FAZ,BAZ,S)
C
C *** SOLUTION OF THE GEODETIC INVERSE PROBLEM AFTER T.VINCENTY.
C *** MODIFIED RAINSFORD WITH HELMERT ELLIPTICAL TERMS.
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL.
C *** STANDPOINT/FOREPOINT MUST NOT BE THE GEOGRAPHIC POLE .
C *** KODE = +1 FOR U.S. DATUM, KODE = -1 FOR EUROPEAN.
                        
c     I don't know about European and U.S. Datum so use the following convention
c     input:
c     E1, E2        Geodetic Longitudes E positive W negative, in radians
c     GLAT1, GLAT2  Geodetic Latitudes, N postive, S negative, in radians
c     
c     output:
c     FAZ           forward (1 to 2) azimuth from 1 to 2, radians clockwise from N
c     BAZ           backward (2 to 1) azimuth radians clockwise from N
c     S             distance in meters
c
c     K. Feigl April 1988

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*2 (I-N)
      DATA TOL/0.5D-14/

c     these are for NAD27 (Clarke 1866 Ellipsoid)
      a = 6378206.40d0
      f = 1.0d0/294.97869820d0

c     these are for WGS84
c      a = 6378137.0d0
c      f = 1.0d0/298.257222101d0

      pi = 4.0d0 * datan(1.0d0)
      twopi =2.0d0 * pi

      R = 1.0D0-F                          

c     hardwire the sign convention
      glon1 = e1
      glon2 = e2
c      GLON1 = DBLE(KODE)*E1
c      GLON2 = DBLE(KODE)*E2 

      TU1 = R*DSIN(GLAT1)/DCOS(GLAT1)
      TU2 = R*DSIN(GLAT2)/DCOS(GLAT2)
      CU1 = 1.0D0/DSQRT(TU1*TU1+1.0D0)
      SU1 = CU1*TU1
      CU2 = 1.0D0/DSQRT(TU2*TU2+1.0D0)
      S = CU1*CU2
      BAZ = S*TU2
      FAZ = BAZ*TU1
      X = GLON1-GLON2
  100 SX = DSIN(X)
      CX = DCOS(X)
      TU1 = CU2*SX
      TU2 = SU1*CU2*CX-BAZ
      SY = DSQRT(TU1*TU1+TU2*TU2)
      CY = S*CX+FAZ
      Y = DATAN2(SY,CY)
      SA = S*SX/SY
      C2A = -SA*SA+1.0D0
      CZ = FAZ+FAZ
      IF(C2A .GT. 0.0D0)CZ = -CZ/C2A+CY
      E = CZ*CZ*2.0D0-1.0D0
      C = ((-3.0D0*C2A+4.0D0)*F+4.0D0)*C2A*F/16.0D0
      D = X
      X = ((E*CY*C+CZ)*SY*C+Y)*SA
      X = (1.0D0-C)*X*F+GLON1-GLON2
      IF(DABS(D-X) .GT. TOL)GO TO 100
      FAZ = DATAN2(-TU1,TU2)
      BAZ = DATAN2(CU1*SX,BAZ*CX-SU1*CU2)
      IF(KODE .EQ. +1)GO TO 101
      FAZ = PI-FAZ
      BAZ = PI-BAZ
 101  IF(FAZ .LT. 0.0D0)FAZ = FAZ+TWOPI
      IF(BAZ .LT. 0.0D0)BAZ = BAZ+TWOPI
      X = DSQRT((1.0D0/R/R-1.0D0)*C2A+1.0D0)+1.0D0
      X = (X-2.0D0)/X
      C = 1.0D0-X
      C = (X*X/4.0D0+1.0D0)/C
      D = (0.375D0*X*X-1.0D0)*X
      X = E*CY
      S = 1.0D0-E-E
      S = ((((SY*SY*4.0D0-3.0D0)*S*CZ*D/6.0D0-X)*D/4.0D0+CZ)*SY*D+Y)
     1    *C*A*R
      

      RETURN
      END


