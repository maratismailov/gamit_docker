      subroutine getdat(iyear,imonth,iday)

c     return the run time date

      integer date(5)
      integer iyear,imonth,iday
      real*8 secs(2)
      external idate

       call systime (date,secs)
       iyear =  date(1)
       call fix_y2k(iyear)
       imonth = date(2)
       iday = date(3)
       
c   Sun and HP
c      call idate(iymd(1),iymd(2),iymd(3))
c      iyear = iymd(3)
c      call fix_y2k(iyear)
c      imonth = iymd(1)
c      iday = iymd(2)

c#   DEC
c#      call idate(imonth,iday,iyear)

      return
      end
