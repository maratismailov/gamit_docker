      FUNCTION JULDAY(JMONTH,IDAY,IYEAR,mode)
C
C     Compute a Julian Day given month, day, and year
C     M.E.Ash   Oct 1966
C
C     Input year from 1900 (negative if before 1900); valid from 1601-2000
c     In FONDA, only the time difference is meaningful.  We choose year
c     1900 as the time origin.
cmk   modified to make julian date calcs ok beyond 2000
C
c     mode = 1: from month, day of month, year to Julian day
c     mode = 2: from day of year, year to Julian day
c     mode = 3: from month, day, year to day of year
c
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 JULDAY,JD1900,int40,doy
      integer imonth,iday,iyear,mode,montot,jyear,jmonth 
      logical leapyr


      DIMENSION MONTOT(12)
      DATA MONTOT/0,31,59,90,120,151,181,212,243,273,304,334/
      data JD1900/2415020/
c     data JD1900/0/
      data int40/0/
c
      imonth=jmonth
c
      if (mode.le.0.or.mode.gt.3) then
         print*,'Unknown mode: ',mode
         stop
      endif
C
C          IMONTH=MONTH           (1-12)
C          IDAY  =DAY OF MONTH    (1-31)  or day of year (1-366)
C          IYEAR =YEAR SINCE 1900 (NEGATIVE BEFORE 1900)
C
      IF(IYEAR .GT. 1600) THEN
         JYEAR = IYEAR - 1900
      ELSE
         JYEAR = IYEAR
      ENDIF
      
c     test if the year is between 1600-????
      IF(JYEAR.LT.-400) THEN
         print *,'JULDAY: invalid year = ',jyear
         stop
      ENDIF

      if (mode.eq.1) then
         IF((IMONTH.LE.0).OR.(IMONTH.GT.12)) then
             print *,'JULDAY: invalid month = ',imonth
             stop
         endif
         IF((IDAY.LE.0).OR.(IDAY.GT.31)) then
            print *,'JULDAY: invalid day of month = ',iday
            stop
         endif
      endif

      if (mode.eq.2) then
         IF((IDAY.LE.0).OR.(IDAY.GT.366)) then
            print *,'JULDAY: invalid day of year = ',iday
            stop
         endif
      endif

      if (mode.eq.2) then
c        dummy values for julian day calculations
c        later add DOY value 
         doy=iday
         iday=1
         imonth=1 
      endif

      if (mode.eq.3) then
         julday = int40+(montot(imonth)+iday)
         if (imonth.gt.2.and.leapyr(iyear)) julday = julday+1
         goto 100
      endif

c     use full jd calculations
      julday = (1461*(iyear+4800+(imonth-14)/12))/4+
     .    (367*(imonth-2-12*((imonth-14)/12)))/12-
     .    (3*((iyear+4900+(imonth-14)/12)/100))/4+iday-32075

c     now convert to reference year; implies conversion to mjd
      julday=julday-JD1900   
 
      if (mode.eq.2) then
c         add on doy 
          julday=julday+doy 
      endif

 100  continue
      RETURN
      END
