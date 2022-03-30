      SUBROUTINE MDYJUL(IMONTH,IDAY,IYEAR,JDS)
C
C     Obtain month, day, year for given Julian Day
C     T. Forni May 1968
C
      INTEGER*4 JDS,JD1600,JD,NYR,IC,imonth,iday,iyear,mdn,i
      DIMENSION MDN(13)
      DATA MDN/0,31,59,90,120,151,181,212,243,273,304,334,365/
     1       , JD1600/2305447/
C
C     JD = days since 0 January 1600
      JD= JDS - JD1600
      NYR= JD/365
C     IC = number of centuries since 0 January 1600
   16 IC= NYR/100
C
C        Compute days due to leap years
      IDAY = JD - NYR*365 - (NYR-1)/4 + (NYR+99)/100 - (NYR+399)/400 - 1
      IF (IC.NE.0 ) GOTO 20
      IF (NYR.NE.0 ) GOTO 20
      IDAY = IDAY + 1
   20 IF (IDAY.GT.0 ) GOTO 23
      NYR = NYR - 1
      GOTO 16
C     IYEAR = (0 thru 99) year of the century
   23 IYEAR= NYR - IC * 100
      NYR= IYEAR
      IF (NYR.NE.0 ) GOTO 27
      IF (MOD(IC,4).NE.0 ) GOTO 34
      GOTO 30
   27 IF (MOD(NYR,4).NE.0) GOTO 34
   30 IF (IDAY - 60 ) 34,39,32
   32 IDAY= IDAY - 1
   34 DO 36 I=2,13
      IF (IDAY.LE.MDN(I)) GOTO 40
   36 CONTINUE
      stop ' Stop in MDYJUL'
   39 IMONTH= 2
      IDAY= 29
      GOTO 45
   40 IMONTH= I - 1
      IDAY= IDAY - MDN(IMONTH)
   45 RETURN
      END
