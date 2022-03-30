Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE RADDEG(DIN,DOUT)

      implicit none  

      real*8 radian,dout,temp,din
C
      DIMENSION DOUT(3)
      DATA RADIAN/57.29577951308232D0/
C
C      DAINT(X)=X-DMOD(X,1.D0)
C
      TEMP=RADIAN*DIN
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
      RETURN
      END
