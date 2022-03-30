      Subroutine fix_y2k(iyear)

c     Check for 2 or 4-digit years.  If 2-digits, convert to 4 by guessing
c     the date.   Same as check_y2k except no warning printed.  R. King July 1999

      implicit none

      integer*4 iyear

      if( iyear.le.1900) then 
                  
c        earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
         if( iyear.gt.60.and.iyear.le.99 ) then
            iyear = iyear + 1900
         else
            iyear = iyear + 2000
         endif  

      endif

      return
      end
