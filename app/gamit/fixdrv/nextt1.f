      SUBROUTINE  NEXTT1(IY,IMO,ID,IH,IM,IS,JDELTA,
     *                   JY,JMO,JD,JH,JM,JS)
C***** ROUTINE TO CALCULATE INCREMENTED TIME POINT
C***** INPUT
C*****   IY,IMO,ID,IH,IM,IS : INITIAL TIME POINT
C*****   JDELTA : INCREMENT VALUE (SEC)
C***** OUTPUT
C*****   JY,JMO,JD,JH,JM,JS : INCREMENTED TIME POINT
C***** SUBROUTINE USED...LEAPYR (library)
C
      logical leapyr

      integer*4 id,jd,ih,jh,im,jm,is,js,iy,jy,jdelta
     .        , imo,jmo,jdmax,moday

      DIMENSION  MODAY(12)

      DATA  MODAY/31,28,31,30,31,30,31,31,30,31,30,31/

      JY=IY
      JMO=IMO
      JD=ID
      JH=IH
      JM=IM
      JS=IS
C***** INCREMENT SEC
      JS=JS+JDELTA
   10 CONTINUE
      IF(JS.LE.59)  GOTO 90
      JS=JS-60
      JM=JM+1
C***** INCREMENT MIN
      IF(JM.LE.59)  GO TO 10
      JM=JM-60
      JH=JH+1
C***** INCREMENT HOUR
      IF(JH.LE.23)  GO TO 10
      JH=JH-24
      JD=JD+1
C***** INCREMENT DAY
      JDMAX=MODAY(JMO)
      IF( JMO.EQ.2. and.leapyr(jy) ) JDMAX=JDMAX+1
      IF(JD.LE.JDMAX)  GO TO 10
      JD=JD-JDMAX
      JMO=JMO+1
C***** INCREMENT MONTH
      IF(JMO.LE.12)  GO TO 10
      JMO=JMO-12
C***** INCREMENT YEAR
      JY=JY+1
      GO TO 10
C*****
   90 RETURN
      END
