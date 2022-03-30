      subroutine conv_date(yr,mo,day,hr,min,sec,date)

c Purpose:  to take a year,mo,day,hr,min and create the strange character
c           string that is used in sinex files for representing this
c
c P Tregoning
c 8th August, 1995

      implicit none

      integer yr,yr2,mo,day,hr,min
      character*12 date
      real*8 sec,doy,sod
      
      yr2 = mod(yr,100)

c  convert to day of year and seconds of day
      call daynum(mo,day,yr,doy)
      sod = sec + min*60.d0 + hr*3600.d0

c  now form up the date string
      write(date,'(i2,a1,i3.3,a1,i5.5)')yr,':',int(doy),':'
     .               ,int(sod)

      return
      end
