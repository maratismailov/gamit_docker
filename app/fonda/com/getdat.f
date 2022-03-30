      subroutine getdat(iyear,imonth,iday) 

c     return the run time date
c     Sun version

      integer i3(3)   
      integer iyear,imonth,iday

      call idate(i3)

      iyear = i3(3)
      imonth = i3(2)
      iday = i3(1)
      return
      end
