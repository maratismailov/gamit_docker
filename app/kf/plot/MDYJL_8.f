CTITLE    ................................................................
 
C@MDYJL_8
      SUBROUTINE MDYJL_8(IMONTH,IDAY,IYEAR,ITIME,FJD)
C
c mod tah 850617 modified to use real*8 variables for the SOLVK
c     program.
c
C      *****************************************************************
C      * MDYJL CONVERTS JULIAN DATE @ MIDNIGHT TO MONTH,DAY, AND YEAR. *
C      * COPIED FROM PEP PLUS AN ADDITIONAL .5 DAY TO CONVERT FROM     *
C      * JULIAN DATE AT MIDNIGHT TO THE PEP JULIAN DAY NUMBER          *
C      *                                                               *
C      *****************************************************************
C
      real*8 FJD, xjd
 
      integer*4 MDN(13)
      integer*4 i, imonth, iday, iyear, itime, ic, nyr, inyr
 
c
      DATA MDN/0,31,59,90,120,151,181,212,243,273,304,334,365/
C
C**** XJD = DAYS SINCE 0 JANUARY, 1600
      XJD = FJD - 2305447.0D0 + 0.5D0
      NYR = XJD/365.
C**** IC = NUMBER OF CENTURIES SINCE 0 JANUARY, 1600
   16 IC = NYR/100
C     DAYS DUE TO LEAP YEARS
      INYR = XJD - NYR*365.0
      IDAY = INYR - (NYR-1)/4 + (NYR + 99)/100 - (NYR + 399)/400 - 1
      IF(IC .NE.0) GO TO 20
      IF(NYR.NE.0) GO TO 20
      IDAY = IDAY + 1
   20 IF(IDAY .GT. 0) GO TO 23
      NYR = NYR - 1
      GO TO 16
C**** IYEAR (O THRU 99) YEAR OF THE CENTURY
   23 IYEAR = NYR - IC * 100
      ITIME = IC - 3
      NYR = IYEAR
      IF(NYR .NE. 0) GO TO 27
      IF(MOD(IC,4) .NE. 0) GO TO 34
      GO TO 30
   27 IF(MOD(NYR,4) .NE. 0) GO TO 34
   30 IF(IDAY - 60) 34,39,32
   32 IDAY = IDAY - 1
   34 DO 36 I=2,13
          IF(IDAY .LE. MDN(I)) GO TO 40
   36 CONTINUE
   39 IMONTH = 2
      IDAY = 29
      GO TO 45
   40 IMONTH = I - 1
      IDAY = IDAY - MDN(IMONTH)
   45 RETURN
      END
 
