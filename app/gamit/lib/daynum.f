      SUBROUTINE DAYNUM(IMONTH,IDAY,IYEAR,DOY)
C
C     DETERMINES THE DAY OF YEAR FROM THE CALENDAR MONTH AND DAY
C     RICK ABBOT - NOVEMBER 1984

      implicit none

      integer*4 month,ly,iday,isub,idoy,iy4,imonth,iyear,iy
      real*8 doy

      DIMENSION MONTH(13,2)

      DATA MONTH /0,31,60,91,121,152,182,213,244,274,305,335,366,
     $            0,31,59,90,120,151,181,212,243,273,304,334,365/

      ISUB(IY) = MIN0( MOD(IY,4),1) + 1
      IY4 = IYEAR
      LY = ISUB(IY4)
      IDOY = MONTH(IMONTH,LY) + IDAY
      DOY=DBLE(IDOY)
C
      RETURN
      END
