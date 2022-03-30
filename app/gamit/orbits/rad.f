      SUBROUTINE RAD (IDEG,MIN,SEC,ANGLE)
C     ******************************************************************
C     **                                                              **
C     **      CONVERTS ANGLE IN DEGREES,MINUTES,SECONDS TO RADIANS    **
C     **                                                              **
C     ******************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      integer*4 ideg,min
      real*8 pi,angle
      pi=4.d0*datan(1.d0)
C
      ANGLE=IABS(IDEG)+IABS(MIN)/60.D0+DABS(SEC)/3600.D0
      IF((IDEG.LT.0).OR.(MIN.LT.0).OR.(SEC.LT.0.D0)) ANGLE=-ANGLE
      ANGLE=ANGLE*PI/180.D0
C
      RETURN
      END
