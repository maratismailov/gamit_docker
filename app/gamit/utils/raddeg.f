Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE RADDEG(DIN,DOUT)

c NOTES: routine from SOLVE
c modified by PT 950810: - implicit none added
c                        - will return the correct sign for dout

      IMPLICIT none
      real*8 din,dout(3),radian,temp
      integer sgn
C
      DATA RADIAN/57.29577951308232D0/
C
C      DAINT(X)=X-DMOD(X,1.D0)
C
c PT950810: save the original sign of the coordinate
      sgn = sign(1.d0,din)
      TEMP=dabs(RADIAN*DIN)
C
CD     WRITE(6,2345) DIN,RADIAN,TEMP
CD2345 FORMAT(1X,D30.15)
C
C     TEMP=TEMP+1.D-9 !SO THAT INTEGER DEGREES WILL
C     FORMAT RIGHT
      DOUT(1)=DINT(TEMP)
      TEMP=60.D0*(TEMP-DOUT(1))
      DOUT(2)=DINT(TEMP)
      DOUT(3)=60.D0*(TEMP-DOUT(2))
C
c PT950810: assign the correct sign to the output degrees
      dout(1) = sign(dout(1),dble(sgn))
      RETURN
      END
