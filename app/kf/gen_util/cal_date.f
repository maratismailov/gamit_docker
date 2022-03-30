ctitle
 
      subroutine cal_date(epoch,date)

      implicit none 
 
c
c     subroutine to convert epoch to the the nearest calender date
c     and hour.
c
c Variables
c ---------
c epoch -- julian day number plus fractional part of day
c date  -- calender date yr,mon,day, hour, min
c
      real*8 epoch
 
c
      integer*4 date(5)
 
c
c Local variables
c ---------------
C sectag - Seconds tag (if greater then 30 seconds we set to next
c          minutes.
c new_epoch - epoch modified so that to nearest miniute
c
      real*8 sectag, new_epoch
 
c
*     Convert Julian date to calender date ans seconds
      call jd_to_ymdhms( epoch, date, sectag )
c
c.... Truncate to nearest minute
      if( sectag.ge.30 ) then
         date(5) = date(5) + 1
c
c....    Now recompute julian day
         sectag = 0.01d0
         call ymdhms_to_jd( date, sectag, new_epoch )

*        Now convert back
         call jd_to_ymdhms( new_epoch, date, sectag )
      end if
c
c**** Thats all
      return
      end
 
c........................................................................
