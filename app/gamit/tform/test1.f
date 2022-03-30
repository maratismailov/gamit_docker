      double precision alat,along,rad,x(3)

      read(5,*) alat,along,rad
      call sphxyz(alat,along,rad,x)
      print*,x(1),x(2),x(3)
      end
      SUBROUTINE SPHXYZ(ALAT,ALONG,RAD,X)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(3)
C
      PI= 4.D0*DATAN(1.D0)
C
      ALAT = ALAT*PI/180.D0
      ALONG = ALONG*PI/180.D0
      X(1)= RAD*DCOS(ALAT)*DCOS(ALONG)
      X(2)= RAD*DCOS(ALAT)*DSIN(ALONG)
      X(3)= RAD*DSIN(ALAT)
C
      RETURN
      END
