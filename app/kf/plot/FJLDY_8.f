CTITLE    ................................................................
 
      real*8 FUNCTION FJLDY_8(IMONTH,IDAY,IYEAR)
 
C
C**** FUNCTION TO CALCULATE THE JULIAN DATE AT MIDNIGHT
C**** THIS ROUTINE FAILS AT 2100 A.D.
C
      real*8 STJD, FYR, YRDY
 
      integer*4 imonth, iday, iyear, montot(12), iyr1, iyr
 
      DATA MONTOT/0,31,59,90,120,151,181,212,243,273,304,334/
C
C**** IDAY = DAY OF MONTH (1-31)
C**** IMONTH = MONTH (1-12)
C**** IYEAR = YEAR SINCE 1900 (NEGATIVE BEFORE 1900)
C
***** See if year is greater than 1000
      if( iyear.gt.1000 ) iyear = iyear - 1900
 
      IYR1 = 0
      STJD = 2415020.0D0
      FYR = IYEAR
      YRDY = 365.0
      IYR = IYEAR/4
      IF(IYEAR) 3,21,7
    3 IYR1 = IYEAR/100
      IF(IYEAR.NE.IYR1*100) GO TO 7
      IF(IMONTH.GT.2) IYR1=IYR1+1
    7 IF(IYEAR.NE.IYR*4) GO TO 21
      IF(IYR) 11,21,15
   11 IF(IMONTH.LE.2) GO TO 21
      IYR = IYR +1
      GO TO 21
   15 IF(IMONTH.GT.2) GO TO 21
      IYR=IYR-1
   21 FJLDY_8 = (MONTOT(IMONTH)+IDAY+IYR-IYR1)
      fjldy_8=fjldy_8 + STJD + 365.*FYR
      fjldy_8 = fjldy_8 - 0.5D0
      RETURN
      END
 
