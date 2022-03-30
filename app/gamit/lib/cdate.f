      subroutine cdate (iyear,idoy,imonth,iday)
c     calculate calendar date from day of year

      integer  iyear,idoy,imonth,iday,imonths(13,2),ly,isub,i,iy


      data imonths/0,31,60,91,121,152,182,213,244,274,305,335,366,
     .             0,31,59,90,120,151,181,212,243,273,304,334,365/

      isub(iy) = min0(mod(iy,4),1) + 1

      ly = isub(iyear)

c     find month
      do 10 i=1,12
         if (idoy .le. imonths(i+1,ly)) goto 20
 10   continue
 20   iday = idoy - imonths(i,ly)
      imonth = i

      return
      end
